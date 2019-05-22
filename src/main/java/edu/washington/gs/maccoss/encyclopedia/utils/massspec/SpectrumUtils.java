package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class SpectrumUtils {
	public static Spectrum getSimpleSpectrum(String name, double precursorMz, float rtInSec, double[] mass, float[] intensity) {
		final float tic=General.sum(intensity);
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
				return precursorMz;
			}
			
			@Override
			public double[] getMassArray() {
				return mass;
			}
			
			@Override
			public float[] getIntensityArray() {
				return intensity;
			}
		};
	}
	public static Spectrum mergeSpectra(ArrayList<Spectrum> spectra, MassTolerance tolerance) {
		if (spectra.size()>50) {
			return binnedMergeSpectra(spectra, 0.1f);
		} else {
			return accurateMergeSpectra(spectra, tolerance);
		}
	}
	public static Spectrum binnedMergeSpectra(ArrayList<Spectrum> spectra, double binWidth) {
		double maxMz=0.0;
		for (Spectrum spectrum : spectra) {
			double mz=spectrum.getMassArray()[spectrum.getMassArray().length-1];
			if (maxMz<mz) maxMz=mz;
		}
		float[] bins=new float[(int)Math.ceil(maxMz/binWidth)];
		if (bins.length==0) return  new PrecursorScan("Combined", 0, 0.0f, null, new double[0], new float[0], 0.0f);

		float totalIIT=0.0f;
		float minRT=Float.MAX_VALUE;
		float tic=0f;
		for (Spectrum spectrum : spectra) {
			if (spectrum.getScanStartTime()<minRT) minRT=spectrum.getScanStartTime();
			if (spectrum instanceof AcquiredSpectrum) {
				float iit=((AcquiredSpectrum)spectrum).getIonInjectionTime();
				if (iit>0) {
					totalIIT+=iit;
				}
			}
			
			double[] mz=spectrum.getMassArray();
			float[] intens=spectrum.getIntensityArray();
			
			for (int i=0; i<mz.length; i++) {
				int index=(int)Math.round(mz[i]/binWidth);
				if (index<0) index=0;
				if (index>=bins.length) index=bins.length-1;
				bins[index]+=intens[i];
			}
			tic += spectrum.getTIC();
		}

		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		for (int i=0; i<bins.length; i++) {
			if (bins[i]>0.0f) {
				masses.add(i*binWidth);
				intensities.add(bins[i]);
			}
		}
		
		return new PrecursorScan("Combined", 0, minRT, totalIIT, masses.toArray(), intensities.toArray(), tic);
	}
	public static Spectrum accurateMergeSpectra(ArrayList<Spectrum> spectra, MassTolerance tolerance) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();

		float totalIIT=0.0f;
		float minRT=Float.MAX_VALUE;
		float tic=0f;
		for (Spectrum spectrum : spectra) {
			if (spectrum.getScanStartTime()<minRT) minRT=spectrum.getScanStartTime();
			if (spectrum instanceof AcquiredSpectrum) {
				float iit=((AcquiredSpectrum)spectrum).getIonInjectionTime();
				if (iit>0) {
					totalIIT+=iit;
				}
			}
			
			double[] mz=spectrum.getMassArray();
			float[] intens=spectrum.getIntensityArray();
			
			for (int i=0; i<mz.length; i++) {
				int index=getIndex(masses, mz[i], tolerance);
				if (index<0) {
					int insertionPoint=-(index+1);
					masses.insert(insertionPoint, mz[i]);
					intensities.insert(insertionPoint, intens[i]);
				} else {
					intensities.setQuick(index, intensities.getQuick(index)+intens[i]);
				}
			}
			tic += spectrum.getTIC();
		}
		return new PrecursorScan("Combined", 0, minRT, totalIIT, masses.toArray(), intensities.toArray(), tic);
	}

	public static int getIndex(TDoubleArrayList peaks, double target, MassTolerance tolerance) {
		if (peaks.size()==0) return -1;
		
		int value=peaks.binarySearch(target);
		// exact match (not likely)
		if (value>=0) return value;
		
		int insertionPoint=-(value+1);
		
		if (insertionPoint>0) {
			// look below
			if (tolerance.compareTo(peaks.get(insertionPoint-1), target)==0) {
				return insertionPoint-1;
			}
		}
		if (insertionPoint<peaks.size()) {
			// look up
			if (tolerance.compareTo(peaks.get(insertionPoint), target)==0) {
				return insertionPoint;
			}
		}
		
		return value;
	}
	
	public static String toDTAString(Spectrum s) {
		StringBuilder sb=new StringBuilder();
		sb.append(s.getSpectrumName());
		sb.append('\n');
		
		sb.append(s.getPrecursorMZ());
		sb.append('\t');
		sb.append(s.getScanStartTime());
		sb.append('\t');
		sb.append(s.getTIC());
		sb.append('\n');
		
		for (int i=0; i<s.getMassArray().length; i++) {
			sb.append(s.getMassArray()[i]);
			sb.append('\t');
			sb.append(s.getIntensityArray()[i]);
			sb.append('\n');
		}
		return sb.toString();
	}
}
