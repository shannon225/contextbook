package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;

public class MedianChromatogramData {
	private final ArrayList<float[]> normalizedChromatograms;
	private final IntRange indices;
	private final float[] medianChromatogram;
	public MedianChromatogramData(ArrayList<float[]> normalizedChromatograms, IntRange indices, float[] medianChromatogram) {
		this.normalizedChromatograms = normalizedChromatograms;
		this.indices = indices;
		this.medianChromatogram = medianChromatogram;
	}
	
	public ArrayList<float[]> getNormalizedChromatograms() {
		return normalizedChromatograms;
	}
	
	public IntRange getIndices() {
		return indices;
	}
	
	public float[] getMedianChromatogram() {
		return medianChromatogram;
	}
}
