package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class SpectrumUtils {
	public static Spectrum mergeSpectra(List<? extends Spectrum> spectra, MassTolerance tolerance) {
		if (spectra.size()>50) {
			return binnedMergeSpectra(spectra, 0.1f);
		} else {
			return accurateMergeSpectra(spectra, tolerance);
		}
	}
	public static Spectrum binnedMergeSpectra(List<? extends Spectrum> spectra, double binWidth) {
		double maxMz=0.0;
		for (Spectrum spectrum : spectra) {
			double mz=spectrum.getMassArray()[spectrum.getMassArray().length-1];
			if (maxMz<mz) maxMz=mz;
		}
		float[] bins=new float[(int)Math.ceil(maxMz/binWidth)];
		if (bins.length==0) return  new PrecursorScan("Combined", 0, 0.0f, 0, 0.0, Double.MAX_VALUE, null, new double[0], new float[0], Optional.empty(), 0.0f);

		float totalIIT=0.0f;
		float minRT=Float.MAX_VALUE;
		
		boolean gotAcquiredData=false;
		float tic=0f;
		int minFraction=Integer.MAX_VALUE;
		double isolationWindowLower=Double.MAX_VALUE;
		double isolationWindowUpper=0.0;
		
		for (Spectrum spectrum : spectra) {
			if (spectrum.getScanStartTime()<minRT) minRT=spectrum.getScanStartTime();
			if (spectrum instanceof AcquiredSpectrum) {
				gotAcquiredData=true;
				AcquiredSpectrum acquiredSpectrum = (AcquiredSpectrum)spectrum;
				float iit=acquiredSpectrum.getIonInjectionTime();
				if (iit>0) {
					totalIIT+=iit;
				}

				if (acquiredSpectrum.getFraction()<minFraction) minFraction=acquiredSpectrum.getFraction();
				if (acquiredSpectrum.getIsolationWindowLower()<isolationWindowLower) isolationWindowLower=acquiredSpectrum.getIsolationWindowLower();
				if (acquiredSpectrum.getIsolationWindowUpper()<isolationWindowUpper) isolationWindowUpper=acquiredSpectrum.getIsolationWindowUpper();
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
		
		if (!gotAcquiredData) {
			minFraction=0;
			isolationWindowLower=0.0;
			isolationWindowUpper=Double.MAX_VALUE;
		}

		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		for (int i=0; i<bins.length; i++) {
			if (bins[i]>0.0f) {
				masses.add(i*binWidth);
				intensities.add(bins[i]);
			}
		}
		
		return new PrecursorScan("Combined", 0, minRT, minFraction, isolationWindowLower, isolationWindowUpper, totalIIT, masses.toArray(), intensities.toArray(), Optional.empty(), tic);
	}
	public static Spectrum accurateMergeSpectra(List<? extends Spectrum> spectra, MassTolerance tolerance) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();

		float totalIIT=0.0f;
		float minRT=Float.MAX_VALUE;
		
		boolean gotAcquiredData=false;
		float tic=0f;
		int minFraction=Integer.MAX_VALUE;
		double isolationWindowLower=Double.MAX_VALUE;
		double isolationWindowUpper=0.0;
		
		for (Spectrum spectrum : spectra) {
			if (spectrum.getScanStartTime()<minRT) minRT=spectrum.getScanStartTime();
			if (spectrum instanceof AcquiredSpectrum) {
				gotAcquiredData=true;
				AcquiredSpectrum acquiredSpectrum = (AcquiredSpectrum)spectrum;
				float iit=acquiredSpectrum.getIonInjectionTime();
				if (iit>0) {
					totalIIT+=iit;
				}

				if (acquiredSpectrum.getFraction()<minFraction) minFraction=acquiredSpectrum.getFraction();
				if (acquiredSpectrum.getIsolationWindowLower()<isolationWindowLower) isolationWindowLower=acquiredSpectrum.getIsolationWindowLower();
				if (acquiredSpectrum.getIsolationWindowUpper()<isolationWindowUpper) isolationWindowUpper=acquiredSpectrum.getIsolationWindowUpper();
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
		
		if (!gotAcquiredData) {
			minFraction=0;
			isolationWindowLower=0.0;
			isolationWindowUpper=Double.MAX_VALUE;
		}
		
		return new PrecursorScan("Combined", 0, minRT, minFraction, isolationWindowLower, isolationWindowUpper, totalIIT, masses.toArray(), intensities.toArray(), Optional.empty(), tic);
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
