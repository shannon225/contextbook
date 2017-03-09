package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;

public class XCorrLibraryEntry extends LibraryEntry {
	private final boolean isDecoy;
	private final SearchParameters params;
	private final SparseXCorrSpectrum spectrum;
	private SparseXCorrCalculator xcorrSpectrum=null;

	public XCorrLibraryEntry(boolean isDecoy, String source, HashSet<String> accessions, byte precursorCharge, String peptideModSeq,
			SparseXCorrSpectrum spectrum, SearchParameters params) {
		super(source, accessions, spectrum.getPrecursorMZ(), precursorCharge, peptideModSeq, 1, 0.0f, // (float)SSRCalc.getHydrophobicity(peptideModSeq)
				0.0f, spectrum.getMassArray(), spectrum.getIntensityArray());
		this.isDecoy=isDecoy;
		this.spectrum=spectrum;
		this.params=params;
	}
	
	public static XCorrLibraryEntry generateEntry(boolean isDecoy, String source, HashSet<String> accessions, byte precursorCharge, String peptideModSeq, SearchParameters params) {
		SparseXCorrSpectrum spectrum=SparseXCorrCalculator.getTheoreticalSpectrum(peptideModSeq, precursorCharge, params);
		return new XCorrLibraryEntry(isDecoy, source, accessions, precursorCharge, peptideModSeq, spectrum, params);
	}
	
	public void init() {
		if (xcorrSpectrum==null) {
			this.xcorrSpectrum=new SparseXCorrCalculator(spectrum, params);
		}
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
