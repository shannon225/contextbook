package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;

public class XCorrLibraryEntry extends PecanLibraryEntry {
	private final SparseXCorrCalculator xcorrSpectrum;

	public XCorrLibraryEntry(FastaPeptideEntry entry, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, boolean isDecoy, SparseXCorrCalculator xcorrSpectrum) {
		super(entry, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, isDecoy, 0.0f);
		this.xcorrSpectrum=xcorrSpectrum;
	}

	public float score(SparseXCorrSpectrum spectrum) {
		return xcorrSpectrum.score(spectrum);
	}
}
