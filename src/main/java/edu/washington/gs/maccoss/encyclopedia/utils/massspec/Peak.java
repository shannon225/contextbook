package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Collection;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.PointInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class Peak implements PointInterface {
	public final double mass;
	public final float intensity;
	public final Float ionMobility;

	public Peak(double mass, float intensity) {
		this(mass, intensity, null);
	}
	public Peak(double mass, float intensity, Float ionMobility) {
		this.mass = mass;
		this.intensity = intensity;
		this.ionMobility=ionMobility;
	}
	
	@Override
	public String toString() {
		return "("+mass+","+intensity+")";
	}

	@Override
	public double getX() {
		return mass;
	}

	@Override
	public double getY() {
		return intensity;
	}
	
	public float getIntensity() {
		return intensity;
	}
	public double getMass() {
		return mass;
	}
	public Float getIonMobility() {
		return ionMobility;
	}
	
	/**
	 * doesn't compare on Y (intensity)
	 */
	@Override
	public int compareTo(PointInterface o) {
		if (o==null) return 1;
		if (mass>o.getX()) return 1;
		if (mass<o.getX()) return -1;
		return 0;
	}
	
	public static Triplet<double[], float[], Optional<float[]>> toArrays(Collection<Peak> peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		TFloatArrayList ionMobilities=new TFloatArrayList();
		for (Peak peak : peaks) {
			masses.add(peak.mass);
			intensities.add(peak.intensity);
			if (peak.ionMobility!=null) {
				ionMobilities.add(peak.ionMobility);
			}
		}
		
		return new Triplet<double[], float[], Optional<float[]>>(masses.toArray(), intensities.toArray(), ionMobilities.size()==0?Optional.empty():Optional.of(ionMobilities.toArray()));
	}
	
	public static Triplet<double[], float[], Optional<float[]>> toArrays(Peak[] peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		TFloatArrayList ionMobilities=new TFloatArrayList();
		for (Peak peak : peaks) {
			masses.add(peak.mass);
			intensities.add(peak.intensity);
			if (peak.ionMobility!=null) {
				ionMobilities.add(peak.ionMobility);
			}
		}
		return new Triplet<double[], float[], Optional<float[]>>(masses.toArray(), intensities.toArray(), ionMobilities.size()==0?Optional.empty():Optional.of(ionMobilities.toArray()));
	}
	
	public static Spectrum toSpectrum(Peak[] peaks, final String name, final float rtInSec, final double precursorMZ) {
		final double[] masses=new double[peaks.length];
		final float[] intensities=new float[peaks.length];
		float[] ims=new float[peaks.length];
		boolean anyIMS=false;
		for (int i = 0; i < peaks.length; i++) {
			masses[i]=peaks[i].mass;
			intensities[i]=peaks[i].intensity;
			if (peaks[i].ionMobility!=null) {
				ims[i]=peaks[i].ionMobility;
				anyIMS=true;
			}
		}
		final float[] finalIMS=anyIMS?ims:null;
		final float tic=General.sum(intensities);
		
		return new Spectrum() {
			
			@Override
			public float getTIC() {
				return tic;
			}
			
			@Override
			public String getSpectrumName() {
				return name;
			}
			
			@Override
			public float getScanStartTime() {
				return rtInSec;
			}
			
			@Override
			public double getPrecursorMZ() {
				return precursorMZ;
			}
			
			@Override
			public double[] getMassArray() {
				return masses;
			}
			
			@Override
			public Optional<float[]> getIonMobilityArray() {
				return Optional.ofNullable(finalIMS);
			}
			
			@Override
			public float[] getIntensityArray() {
				return intensities;
			}
		};
	}
}
