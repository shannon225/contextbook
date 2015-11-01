package edu.washington.gs.maccoss.encyclopedia.datastructures;

public class PrecursorScan {
	private final String spectrumName;
	private final int spectrumIndex;
	private final float scanStartTime;
	private final double[] massArray;
	private final float[] intensityArray;

	public PrecursorScan(String spectrumName, int spectrumIndex, float scanStartTime, double[] massArray, float[] intensityArray) {
		this.spectrumName=spectrumName;
		this.spectrumIndex=spectrumIndex;
		this.scanStartTime=scanStartTime;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
	}

	public String getSpectrumName() {
		return spectrumName;
	}

	public int getSpectrumIndex() {
		return spectrumIndex;
	}

	public float getScanStartTime() {
		return scanStartTime;
	}

	public double[] getMassArray() {
		return massArray;
	}

	public float[] getIntensityArray() {
		return intensityArray;
	}

}
