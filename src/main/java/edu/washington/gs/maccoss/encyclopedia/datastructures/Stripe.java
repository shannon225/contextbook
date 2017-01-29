package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

//@Immutable
public class Stripe implements Comparable<Stripe>, Spectrum {
	private final String spectrumName;
	private final String precursorName;
	private final int spectrumIndex;
	private final float scanStartTime;
	private final float isolationWindowLower;
	private final float isolationWindowUpper;
	private final double[] massArray;
	private final float[] intensityArray;
	private final float intensityMagnitude;
	private final float tic;

	public Stripe(String spectrumName, String precursorName, int spectrumIndex, float scanStartTime, float isolationWindowLower, float isolationWindowUpper, double[] massArray, float[] intensityArray) {
		this.spectrumName=spectrumName;
		this.precursorName=precursorName;
		this.spectrumIndex=spectrumIndex;
		this.scanStartTime=scanStartTime;
		this.isolationWindowLower=isolationWindowLower;
		this.isolationWindowUpper=isolationWindowUpper;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
		
		float thisTic=0.0f;
		float magnitude=0.0f;
		for (float f : intensityArray) {
			thisTic+=f;
			magnitude+=f*f;
		}
		intensityMagnitude=(float)Math.sqrt(magnitude);
		tic=thisTic;
	}
	
	public Stripe sqrt() {
		return new Stripe(spectrumName, precursorName, spectrumIndex, scanStartTime, isolationWindowLower, isolationWindowUpper, massArray, General.protectedSqrt(intensityArray));
	}
	
	@Override
	public int compareTo(Stripe o) {
		if (o==null) return 1;
		int c=Float.compare(scanStartTime, o.scanStartTime);
		if (c!=0) return c;
		c=Integer.compare(spectrumIndex, o.spectrumIndex);
		if (c!=0) return c;
		c=Float.compare(isolationWindowLower, o.isolationWindowLower);
		if (c!=0) return c;
		c=Float.compare(isolationWindowUpper, o.isolationWindowUpper);
		if (c!=0) return c;
		return spectrumName.compareTo(o.spectrumName);
	}
	
	public float getIntensityMagnitude() {
		return intensityMagnitude;
	}
	
	public float getTIC() {
		return tic;
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
	
	public float getIsolationWindowCenter() {
		return (isolationWindowLower+isolationWindowUpper)/2.0f;
	}
	
	public double getPrecursorMZ() {
		return getIsolationWindowCenter();
	}

	public double[] getMassArray() {
		return massArray;
	}

	public float[] getIntensityArray() {
		return intensityArray;
	}
	
	public static ArrayList<Spectrum> downcastStripeToSpectrum(ArrayList<Stripe> stripes) {
		ArrayList<Spectrum> spectra=new ArrayList<Spectrum>();
		for (Spectrum spectrum : stripes) {
			spectra.add(spectrum);
		}
		return spectra;
	}
}
