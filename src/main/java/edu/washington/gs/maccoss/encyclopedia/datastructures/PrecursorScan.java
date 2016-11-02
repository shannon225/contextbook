package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

//@Immutable
public class PrecursorScan implements Spectrum, Comparable<PrecursorScan> {
	private final String spectrumName;
	private final int spectrumIndex;
	private final float scanStartTime;
	private final double[] massArray;
	private final float[] intensityArray;
	private final float tic;

	public PrecursorScan(String spectrumName, int spectrumIndex, float scanStartTime, double[] massArray, float[] intensityArray) {
		this.spectrumName=spectrumName;
		this.spectrumIndex=spectrumIndex;
		this.scanStartTime=scanStartTime;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
		this.tic=General.sum(intensityArray);
	}
	
	@Override
	public float getTIC() {
		return tic;
	}
	
	@Override
	public double getPrecursorMZ() {
		return -1.0;
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

	@Override
	public int compareTo(PrecursorScan o) {
		if (o==null) return 1;
		int c=Float.compare(scanStartTime, o.scanStartTime);
		if (c!=0) return c;

		c=Integer.compare(spectrumIndex, o.spectrumIndex);
		if (c!=0) return c;

		c=spectrumName.compareTo(o.spectrumName);
		return c;
	}
}
