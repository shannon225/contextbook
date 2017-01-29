package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;

public class XCorrLibraryEntry extends PecanLibraryEntry {
	private final SparseXCorrCalculator xcorrSpectrum;

	public XCorrLibraryEntry(FastaPeptideEntry entry, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, boolean isDecoy, SearchParameters params) {
		super(entry, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, isDecoy, 0.0f);
		this.xcorrSpectrum=new SparseXCorrCalculator(peptideModSeq, params);
	}

	public float score(SparseXCorrSpectrum spectrum) {
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
