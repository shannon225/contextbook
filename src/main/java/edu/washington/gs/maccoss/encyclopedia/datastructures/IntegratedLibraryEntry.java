package edu.washington.gs.maccoss.encyclopedia.datastructures;


public class IntegratedLibraryEntry extends LibraryEntry {
	private final Range rtRange;
	
	public IntegratedLibraryEntry(int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, Range rtRange) {
		super(spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray);
		this.rtRange=rtRange;
	}
	
	public Range getRtRange() {
		return rtRange;
	}
}
