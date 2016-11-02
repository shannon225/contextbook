package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public interface Spectrum {
	public String getSpectrumName();
	public float getScanStartTime();
	public double getPrecursorMZ();
	public double[] getMassArray();
	public float[] getIntensityArray();
	public float getTIC();
}
