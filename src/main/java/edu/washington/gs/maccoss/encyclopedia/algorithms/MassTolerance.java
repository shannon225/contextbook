package edu.washington.gs.maccoss.encyclopedia.algorithms;

import gnu.trove.list.array.TDoubleArrayList;

import java.util.Arrays;

import com.google.common.base.Optional;

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
		int value=Arrays.binarySearch(peaks, target);
		// exact match (not likely)
		if (value>=0) return Optional.of(peaks[value]);
		
		// insertion point
		value=-(value+1);
		
		if (value>0) {
			// look below
			if (compareTo(peaks[value-1], target)==0) return Optional.of(peaks[value-1]);
		}
		if (value<peaks.length) {
			// look up
			if (compareTo(peaks[value], target)==0) return Optional.of(peaks[value]);
		}
		
		return Optional.absent();
	}
	
	/**
	 * @param peaks -- assumes sorted array of peaks
	 * @param target
	 * @return all matching masses in range
	 */
	public double[] getMatches(double[] peaks, double target) {
		int value=Arrays.binarySearch(peaks, target);
		// exact match (not likely)
		if (value<0) {
			// insertion point
			value=-(value+1);
		}
		
		TDoubleArrayList matches=new TDoubleArrayList();
		// look below
		int index=value;
		while (index>0&&compareTo(peaks[index-1], target)==0) {
			matches.add(peaks[index-1]);
			index--;
		}

		// look up
		index=value;
		while (index<peaks.length&&compareTo(peaks[index], target)==0) {
			matches.add(peaks[index]);
			index++;
		}

		return matches.toArray();
	}
}
