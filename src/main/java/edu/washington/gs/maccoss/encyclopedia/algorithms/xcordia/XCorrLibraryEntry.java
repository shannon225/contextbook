package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;

public class XCorrLibraryEntry extends LibraryEntry {
	private final boolean isDecoy;
	private final SearchParameters params;
	private final SparseXCorrSpectrum spectrum;
	private final FastaPeptideEntry peptide;
	private SparseXCorrCalculator xcorrSpectrum=null;

	public XCorrLibraryEntry(boolean isDecoy, String peptideModSeq, FastaPeptideEntry peptide, byte precursorCharge, SparseXCorrSpectrum spectrum, SearchParameters params) {
		super(peptide.getFilename(), peptide.getAccessions(), spectrum.getPrecursorMZ(), precursorCharge, peptideModSeq, 1, 0.0f, // (float)SSRCalc.getHydrophobicity(peptideModSeq)
				0.0f, spectrum.getMassArray(), spectrum.getIntensityArray(), params.getAAConstants());
		this.peptide=peptide;
		this.isDecoy=isDecoy;
		this.spectrum=spectrum;
		this.params=params;
	}

	public XCorrLibraryEntry(boolean isDecoy, String peptideModSeq, FastaPeptideEntry peptide, byte precursorCharge, SparseXCorrSpectrum spectrum, float retentionTime, SearchParameters params) {
		super(peptide.getFilename(), peptide.getAccessions(), spectrum.getPrecursorMZ(), precursorCharge, peptideModSeq, 1, retentionTime,
				0.0f, spectrum.getMassArray(), spectrum.getIntensityArray(), params.getAAConstants());
		this.peptide=peptide;
		this.isDecoy=isDecoy;
		this.spectrum=spectrum;
		this.params=params;
	}
	
	public static XCorrLibraryEntry generateEntry(boolean isDecoy, LibraryEntry entry, byte precursorCharge, SearchParameters params) {
		SparseXCorrSpectrum spectrum=SparseXCorrCalculator.normalize(entry, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), false, params);
		FastaPeptideEntry peptide=new FastaPeptideEntry(entry.getSource(), entry.getAccessions(), entry.getPeptideModSeq());
		return new XCorrLibraryEntry(isDecoy, entry.getPeptideModSeq(), peptide, precursorCharge, spectrum, entry.getRetentionTime(), params);
	}
	
	public static XCorrLibraryEntry generateEntry(boolean isDecoy, FastaPeptideEntry peptide, byte precursorCharge, SearchParameters params) {
		Pair<FragmentationModel, SparseXCorrSpectrum> theoreticalSpectrumPair=SparseXCorrCalculator.getTheoreticalSpectrumPair(peptide.getSequence(), precursorCharge, params);
		FragmentationModel model=theoreticalSpectrumPair.x;
		SparseXCorrSpectrum spectrum=theoreticalSpectrumPair.y;
		return new XCorrLibraryEntry(isDecoy, model.getModifiedSequence(), peptide, precursorCharge, spectrum, params);
	}
	
	public void init() {
		if (xcorrSpectrum==null) {
			this.xcorrSpectrum=new SparseXCorrCalculator(spectrum, params);
		}
	}
	
	public FastaPeptideEntry getPeptide() {
		return peptide;
	}
	
	@Override
	public boolean isDecoy() {
		return isDecoy;
	}

	public float score(SparseXCorrSpectrum spectrum) {
		if (xcorrSpectrum==null) init();
		
		return xcorrSpectrum.score(spectrum);
	}
	
	public static ArrayList<LibraryEntry> downcast(ArrayList<XCorrLibraryEntry> entries) {
		ArrayList<LibraryEntry> downcast=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : entries) {
			downcast.add(entry);
		}
		return downcast;
	}
}
