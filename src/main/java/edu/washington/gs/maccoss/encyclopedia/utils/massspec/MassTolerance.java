package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Arrays;
import java.util.Optional;

import gnu.trove.list.array.TIntArrayList;

//@Immutable
public class MassTolerance {
	private final double ppmTolerance;
	private final double percent;

	public MassTolerance(double ppmTolerance) {
		this.ppmTolerance = ppmTolerance;
		this.percent=ppmTolerance/1000000.0;
	}
	
	public double getPpmTolerance() {
		return ppmTolerance;
	}
	public double getTolerance(double m) {
		return ppmTolerance/1000000.0*m;
	}
	
	public boolean equals(double m1, double m2) {
		return compareTo(m1, m2)==0;
	}
	
	/**
	 * if first is less, -1, if second is less 1, otherwise 0
	 * @param m1
	 * @param m2
	 * @return
	 */
	public int compareTo(double m1, double m2) {
		double amuTolerance=Math.max(m1, m2)*percent;
		if (m1+amuTolerance<m2) return -1;
		if (m1-amuTolerance>m2) return 1;
		return 0;
	}
	
	/**
	 * @param peaks -- assumes sorted array of peaks
	 * @param target
	 * @return can return null!
	 */
	public Optional<Double> getMatch(double[] peaks, double target) {
		Optional<Integer> index=getIndex(peaks, target);
		if (index.isPresent()) {
			return Optional.of(peaks[index.get()]);
		} else {
			return Optional.empty();
		}
	}

	public Optional<Integer> getIndex(double[] peaks, double target) {
		int value=Arrays.binarySearch(peaks, target);
		// exact match (not likely)
		if (value>=0) return Optional.of(value);
		
		// insertion point
		value=-(value+1);
		
		if (value>0) {
			// look below
			if (compareTo(peaks[value-1], target)==0) return Optional.of(value-1);
		}
		if (value<peaks.length) {
			// look up
			if (compareTo(peaks[value], target)==0) return Optional.of(value);
		}
		
		return Optional.empty();
	}
	
	/**
	 * @param peaks -- assumes sorted array of peaks
	 * @param target
	 * @return all matching masses in range
	 */
	public int[] getIndicies(double[] peaks, double target) {
		int value=Arrays.binarySearch(peaks, target);
		// exact match (not likely)
		if (value<0) {
			// insertion point
			value=-(value+1);
		}
		
		TIntArrayList matches=new TIntArrayList();
		// look below
		int index=value;
		while (index>0&&compareTo(peaks[index-1], target)==0) {
			matches.add(index-1);
			index--;
		}

		// look up
		index=value;
		while (index<peaks.length&&compareTo(peaks[index], target)==0) {
			matches.add(index);
			index++;
		}

		return matches.toArray();
	}
	
	/**
	 * @param peaks -- assumes sorted array of peaks
	 * @param target
	 * @return all matching masses in range
	 */
	public double[] getMatches(double[] peaks, double target) {
		int[] indicies=getIndicies(peaks, target);
		double[] matches=new double[indicies.length];
		for (int i=0; i<indicies.length; i++) {
			matches[i]=peaks[indicies[i]];
		}
		return matches;
	}
	
	/**
	 * assumes targets and masses are in sorted order (and intensities follows masses)
	 * @param masses
	 * @param intensities
	 * @param targets
	 * @return
	 */
	public float[] getIntegratedIntensities(double[] masses, float[] intensities, double[] targets) {
		float[] tics=new float[targets.length];
		if (targets.length==0||masses.length==0) {
			return tics;
		}
		
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			double targetMass=targets[libraryIndex];
			int compare=compareTo(targetMass, masses[spectrumIndex]);
			if (compare==0) {
				tics[libraryIndex]+=intensities[spectrumIndex];
				//libraryIndex++; // could match multiple acquired peaks to the same library peak
				spectrumIndex++;
			} else if (compare>0) {
				spectrumIndex++;
			} else {
				libraryIndex++;
			}
			if (libraryIndex>=targets.length) break;
			if (spectrumIndex>=masses.length) break;
		}
		
		return tics;
	}
	
	/**
	 * @param peaks -- assumes sorted array of peaks
	 * @param target
	 * @return all matching masses in range
	 */
	public float getIntegratedIntensity(double[] masses, float[] intensities, double target) {
		int[] indicies=getIndicies(masses, target);
		float intensity=0.0f;
		for (int i=0; i<indicies.length; i++) {
			intensity+=intensities[indicies[i]];
		}
		return intensity;
	}
	
	/**
	 * @param peaks -- assumes sorted array of peaks
	 * @param target
	 * @return all matching masses in range
	 */
	public float getMaxIntensity(double[] masses, float[] intensities, double target) {
		int[] indicies=getIndicies(masses, target);
		float intensity=0.0f;
		for (int i=0; i<indicies.length; i++) {
			if(intensity<intensities[indicies[i]]) {
				intensity=intensities[indicies[i]];
			}
		}
		return intensity;
	}
}
