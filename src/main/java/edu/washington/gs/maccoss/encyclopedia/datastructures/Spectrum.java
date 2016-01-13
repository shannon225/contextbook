package edu.washington.gs.maccoss.encyclopedia.datastructures;

public interface Spectrum {
	public String getSpectrumName();
	public float getScanStartTime();
	public double[] getMassArray();
	public float[] getIntensityArray();
}
