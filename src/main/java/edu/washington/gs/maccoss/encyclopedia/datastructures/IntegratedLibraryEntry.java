package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

public class IntegratedLibraryEntry extends LibraryEntry {
	private final Range rtRange;
	
	public IntegratedLibraryEntry(HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, Range rtRange) {
		super(accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray);
		this.rtRange=rtRange;
	}
	
	public Range getRtRange() {
		return rtRange;
	}
}
