package edu.washington.gs.maccoss.encyclopedia.datastructures;

public interface Chromatogram extends Spectrum {
	public Range getRtRange();
	public float[] getCorrelationArray();
	public float[] getMedianChromatogram();
}
