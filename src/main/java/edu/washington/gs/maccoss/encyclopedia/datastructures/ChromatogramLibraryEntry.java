package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

public class ChromatogramLibraryEntry extends LibraryEntry implements Chromatogram {
	private final float[] medianChromatogram;
	private final Range range;
	
	public ChromatogramLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, float[] correlationArray, float[] medianChromatogram, Range range) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray);
		this.medianChromatogram=medianChromatogram;
		this.range=range;
	}
	
	public Range getRtRange() {
		return range;
	}
	
	public float[] getMedianChromatogram() {
		return medianChromatogram;
	}
}
