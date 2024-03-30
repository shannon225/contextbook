package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.AcquiredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumPeakFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SpectrumWithCharge;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

//@Immutable
public class FragmentScan implements Comparable<FragmentScan>, AcquiredSpectrum, SpectrumWithCharge {
	private final String spectrumName;
	private final String precursorName;
	private final int spectrumIndex;
	private final float scanStartTime;
	private final double isolationWindowLower;
	private final double isolationWindowUpper;
	private final double[] massArray;
	private final float[] intensityArray;
	private final Optional<float[]> ionMobilityArray;
	private final float intensityMagnitude;
	private final float ionInjectionTime;
	private final float tic;
	private final byte charge;
	private final int fraction; 

	public FragmentScan(String spectrumName, String precursorName, int spectrumIndex, float scanStartTime, int fraction, Float ionInjectionTime, double isolationWindowLower, double isolationWindowUpper, double[] massArray, float[] intensityArray, float[] ionMobilityArray) {
		this(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray, ionMobilityArray, (byte)0);
	}

	public FragmentScan(String spectrumName, String precursorName, int spectrumIndex, float scanStartTime, int fraction, Float ionInjectionTime, double isolationWindowLower, double isolationWindowUpper, double[] massArray, float[] intensityArray, Optional<float[]> ionMobilityArray) {
		this(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray, ionMobilityArray, (byte)0);
	}

	public FragmentScan(String spectrumName, String precursorName, int spectrumIndex, float scanStartTime, int fraction, Float ionInjectionTime, double isolationWindowLower, double isolationWindowUpper, double[] massArray, float[] intensityArray, float[] ionMobilityArray, byte charge) {
		this(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray, Optional.ofNullable(ionMobilityArray), charge);
	}

	public FragmentScan(String spectrumName, String precursorName, int spectrumIndex, float scanStartTime, int fraction, Float ionInjectionTime, double isolationWindowLower, double isolationWindowUpper, double[] massArray, float[] intensityArray, Optional<float[]> ionMobilityArray, byte charge) {
		this.spectrumName=spectrumName;
		this.precursorName=precursorName;
		this.spectrumIndex=spectrumIndex;
		this.scanStartTime=scanStartTime;
		this.fraction=fraction;
		if (ionInjectionTime==null) ionInjectionTime=-1f;
		this.ionInjectionTime=ionInjectionTime;
		
		this.isolationWindowLower=isolationWindowLower;
		this.isolationWindowUpper=isolationWindowUpper;
		this.massArray=massArray;
		this.intensityArray=intensityArray;
		this.ionMobilityArray=ionMobilityArray;
		this.charge=charge;
		
		float thisTic=0.0f;
		float magnitude=0.0f;
		for (float f : intensityArray) {
			thisTic+=f;
			magnitude+=f*f;
		}
		intensityMagnitude=(float)Math.sqrt(magnitude);
		tic=thisTic;
	}
	
	@Override
	public int getFraction() {
		return fraction;
	}
	
	public FragmentScan shallowClone(int fraction, int spectrumIndex) {
		return new FragmentScan(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massArray, intensityArray, ionMobilityArray.orElse(null), charge);
	}
	
	public FragmentScan sqrt() {
		return new FragmentScan(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, massArray, General.protectedSqrt(intensityArray), ionMobilityArray.orElse(null), charge);
	}
	
	/**
	 * select top N peaks in 100 m/z
	 * @param depth
	 * @return
	 */
	public FragmentScan trimToPeakDepth(int depth) {
		// peaks above 2000 are placed in the last bin
		return SpectrumPeakFilter.filterPeaks(this, 100.0, 20, depth);
	}
	
	public FragmentScan trimMasses(Range r) {
		TFloatArrayList ints=new TFloatArrayList();
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList mobilities=new TFloatArrayList();
		for (int i = 0; i < massArray.length; i++) {
			if (r.contains(massArray[i])) {
				ints.add(intensityArray[i]);
				masses.add(massArray[i]);
				if (ionMobilityArray.isPresent()) {
					mobilities.add(ionMobilityArray.get()[i]);
				}
			}
		}
		float[] mobilitiesArray=ionMobilityArray.isPresent()?mobilities.toArray():null;
		return new FragmentScan(spectrumName, precursorName, spectrumIndex, scanStartTime, fraction, ionInjectionTime, isolationWindowLower, isolationWindowUpper, masses.toArray(), ints.toArray(), mobilitiesArray, charge);
	}
	
	@Override
	public float getIonInjectionTime() {
		return ionInjectionTime;
	}
	
	public Optional<float[]> getIonMobilityArray() {
		return ionMobilityArray;
	}
	
	/**
	 * can return 0 (if charge state is unknown)
	 * @return
	 */
	public byte getPrecursorCharge() {
		return charge;
	}
	
	@Override
	public int compareTo(FragmentScan o) {
		if (o==null) return 1;
		int c=Float.compare(scanStartTime, o.scanStartTime);
		if (c!=0) return c;
		c=Integer.compare(spectrumIndex, o.spectrumIndex);
		if (c!=0) return c;
		c=Double.compare(isolationWindowLower, o.isolationWindowLower);
		if (c!=0) return c;
		c=Double.compare(isolationWindowUpper, o.isolationWindowUpper);
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

	public double getIsolationWindowLower() {
		return isolationWindowLower;
	}

	public double getIsolationWindowUpper() {
		return isolationWindowUpper;
	}
	
	public double getIsolationWindowCenter() {
		return (isolationWindowLower+isolationWindowUpper)/2.0;
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

	/**
	 * @deprecated Instead of using this method, refactor usages of its output to us a
	 *             {@code List<? extends Spectrum>} instead of a {@code List<Spectrum>}
	 *             and pass {@code stripes} to it directly.
	 */
	@Deprecated
	public static ArrayList<Spectrum> downcastStripeToSpectrum(ArrayList<FragmentScan> stripes) {
		ArrayList<Spectrum> spectra=new ArrayList<Spectrum>();
		for (Spectrum spectrum : stripes) {
			spectra.add(spectrum);
		}
		return spectra;
	}
}
