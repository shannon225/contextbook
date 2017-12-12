package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class SpectrumUtils {
	public static Spectrum mergeSpectra(ArrayList<Spectrum> spectra, MassTolerance tolerance) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		
		float minRT=Float.MAX_VALUE;
		float tic=0f;
		for (Spectrum spectrum : spectra) {
			if (spectrum.getScanStartTime()<minRT) minRT=spectrum.getScanStartTime();
			
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
		return new PrecursorScan("Combined", 0, minRT, masses.toArray(), intensities.toArray(), tic);
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
