package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.AcquiredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

//@Immutable
public class PrecursorScan implements AcquiredSpectrum, Comparable<PrecursorScan> {
	private final String spectrumName;
	private final int spectrumIndex;
	private final float scanStartTime;
	private final float ionInjectionTime;
	private final double[] massArray;
	private final float[] intensityArray;
	private final Optional<float[]> ionMobilityArray;
	private final float tic;
	private final int fraction; 
	private final double isolationWindowLower; 
	private final double isolationWindowUpper;

	public PrecursorScan(String spectrumName, int spectrumIndex, float scanStartTime, int fraction, double isolationWindowLower, double isolationWindowUpper, Float ionInjectionTime, double[] massArray, float[] intensityArray, float[] ionMobilityArray) {
		this(spectrumName, spectrumIndex, scanStartTime, fraction, isolationWindowLower, isolationWindowUpper, ionInjectionTime, massArray, intensityArray, Optional.ofNullable(ionMobilityArray), null);
	}

	public PrecursorScan(String spectrumName, int spectrumIndex, float scanStartTime, int fraction, double isolationWindowLower, double isolationWindowUpper, Float ionInjectionTime, double[] massArray, float[] intensityArray, Optional<float[]> ionMobilityArray) {
		this(spectrumName, spectrumIndex, scanStartTime, fraction, isolationWindowLower, isolationWindowUpper, ionInjectionTime, massArray, intensityArray, ionMobilityArray, null);
	}
	public PrecursorScan(String spectrumName, int spectrumIndex, float scanStartTime, int fraction, double isolationWindowLower, double isolationWindowUpper, Float ionInjectionTime, double[] massArray, float[] intensityArray, float[] ionMobilityArray, Float tic) {
		this(spectrumName, spectrumIndex, scanStartTime, fraction, isolationWindowLower, isolationWindowUpper, ionInjectionTime, massArray, intensityArray, Optional.ofNullable(ionMobilityArray), tic);
	}

	public PrecursorScan(String spectrumName, int spectrumIndex, float scanStartTime, int fraction, double isolationWindowLower, double isolationWindowUpper, Float ionInjectionTime, double[] massArray, float[] intensityArray, Optional<float[]> ionMobilityArray, Float tic) {
		this.spectrumName=spectrumName;
		this.spectrumIndex=spectrumIndex;
		this.scanStartTime=scanStartTime;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
		this.ionMobilityArray=ionMobilityArray;
		this.fraction=fraction;
		this.isolationWindowLower=isolationWindowLower;
		this.isolationWindowUpper=isolationWindowUpper;

		if (ionInjectionTime==null) {
			this.ionInjectionTime=-1f;
		} else {
			this.ionInjectionTime=ionInjectionTime;
		}
		
		if (tic==null || !Float.isFinite(tic) || tic <= 0) {
			tic = General.sum(intensityArray);
		}
		this.tic=tic;
	}
	
	public PrecursorScan shallowClone(int fraction, int spectrumIndex) {
		return new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, fraction, isolationWindowLower, isolationWindowUpper, ionInjectionTime, massArray, intensityArray, ionMobilityArray.orElse(null), tic);
	}
	
	/**
	 * clone with potentially narrower isolation window to retarget scan
	 * @param fraction
	 * @param spectrumIndex
	 * @param precursorIsolationWindow
	 * @return
	 */
	public PrecursorScan shallowClone(int fraction, int spectrumIndex, Range precursorIsolationWindow) {
		double lowerBound=Math.max(precursorIsolationWindow.getStart(), isolationWindowLower);
		double UpperBound=Math.min(precursorIsolationWindow.getStop(), isolationWindowUpper);

		return new PrecursorScan(spectrumName, spectrumIndex, scanStartTime, fraction, lowerBound, UpperBound, ionInjectionTime, massArray, intensityArray, ionMobilityArray.orElse(null), tic);
	}
	
	@Override
	public int getFraction() {
		return fraction;
	}
	
	public Optional<float[]> getIonMobilityArray() {
		return ionMobilityArray;
	}
	
	public double getIsolationWindowLower() {
		return isolationWindowLower;
	}
	public double getIsolationWindowUpper() {
		return isolationWindowUpper;
	}
	
	@Override
	public float getTIC() {
		return tic;
	}
	
	public float getIonInjectionTime() {
		return ionInjectionTime;
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

	/**
	 * @deprecated Instead of using this method, refactor usages of its output to use a
	 *             {@code List<? extends Spectrum>} instead of a {@code List<Spectrum>}
	 *             and pass {@code precursors} to it directly.
	 */
	@Deprecated
	public static ArrayList<Spectrum> downcast(ArrayList<PrecursorScan> precursors) {
		ArrayList<Spectrum> spectra=new ArrayList<Spectrum>();
		for (Spectrum spectrum : precursors) {
			spectra.add(spectrum);
		}
		return spectra;
	}
}
