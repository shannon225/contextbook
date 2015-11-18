package edu.washington.gs.maccoss.encyclopedia.datastructures;


//@Immutable
public class Stripe {
	private final String spectrumName;
	private final String precursorName;
	private final int spectrumIndex;
	private final float scanStartTime;
	private final float isolationWindowLower;
	private final float isolationWindowUpper;
	private final double[] massArray;
	private final float[] intensityArray;

	public Stripe(String spectrumName, String precursorName, int spectrumIndex, float scanStartTime, float isolationWindowLower, float isolationWindowUpper, double[] massArray, float[] intensityArray) {
		this.spectrumName=spectrumName;
		this.precursorName=precursorName;
		this.spectrumIndex=spectrumIndex;
		this.scanStartTime=scanStartTime;
		this.isolationWindowLower=isolationWindowLower;
		this.isolationWindowUpper=isolationWindowUpper;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
	}
	
	public Range getRange() {
		return new Range(isolationWindowLower, isolationWindowUpper);
	}

	public String getSpectrumName() {
		return spectrumName;
	}

	public String getPrecursorName() {
		return precursorName;
	}

	public int getSpectrumIndex() {
		return spectrumIndex;
	}

	public float getScanStartTime() {
		return scanStartTime;
	}

	public float getIsolationWindowLower() {
		return isolationWindowLower;
	}

	public float getIsolationWindowUpper() {
		return isolationWindowUpper;
	}

	public double[] getMassArray() {
		return massArray;
	}

	public float[] getIntensityArray() {
		return intensityArray;
	}
}
