package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Arrays;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

//@Immutable
public class MassTolerance implements Comparable<MassTolerance> {
	private static final double UNUSED_TOLERANCE=-1.0;
	private final MassErrorUnitType type; // here we treat Resolution as a cosmetic version of PPM
	private final double amuTolerance;
	private final double ppmTolerance;
	private final double percent;

	public MassTolerance(double ppmTolerance) {
		this(ppmTolerance, MassErrorUnitType.PPM);
	}

	public MassTolerance(double tolerance, MassErrorUnitType type) {
		this.type=type;
		
		if (MassErrorUnitType.PPM==type) {
			this.amuTolerance=UNUSED_TOLERANCE;
			this.ppmTolerance=tolerance;
			this.percent=ppmTolerance/1000000.0;
		} else if (MassErrorUnitType.RESOLUTION==type) {
			this.amuTolerance=UNUSED_TOLERANCE;
			this.ppmTolerance=0.5*1000000/tolerance;
			this.percent=ppmTolerance/1000000.0;
		} else {
			this.amuTolerance=tolerance;
			this.ppmTolerance=UNUSED_TOLERANCE;
			this.percent=ppmTolerance/1000000.0;
		}
	}
	
	@Override
	public int compareTo(MassTolerance o) {
		if (o==null) return 1;
		int c=Double.compare(ppmTolerance, o.ppmTolerance);
		if (c!=0) return c;
		return Double.compare(amuTolerance, o.amuTolerance);
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null||!(obj instanceof MassTolerance)) {
			return false;
		}
		return compareTo((MassTolerance)obj)==0;
	}
	
	@Override
	public int hashCode() {
		return Double.hashCode(amuTolerance)*16807*Double.hashCode(ppmTolerance);
	}
	
	@Override
	public String toString() {
		if (MassErrorUnitType.PPM==type) {
			return ppmTolerance+" "+MassErrorUnitType.toString(type);
		} else if (MassErrorUnitType.RESOLUTION==type) {
			return Math.round(0.5*1000000/ppmTolerance)+" "+MassErrorUnitType.toString(type);
		} else {
			return amuTolerance+" "+MassErrorUnitType.toString(type);
		}
	}
	
	public boolean isRelativeTolerance() {
		if (MassErrorUnitType.PPM==type) {
			return true;
		} else if (MassErrorUnitType.RESOLUTION==type) {
			return true;
		} else {
			return false;
		}
	}
	
	public double getPpmTolerance() {
		return ppmTolerance;
	}
	
	public String getUnits() {
		if (amuTolerance==UNUSED_TOLERANCE) {
			return MassErrorUnitType.toString(MassErrorUnitType.PPM);
		} else {
			return MassErrorUnitType.toString(MassErrorUnitType.AMU);
		}
	}
	
	public double getTolerance(double mass) {
		if (amuTolerance==UNUSED_TOLERANCE) {
			return percent*mass;
		} else {
			return amuTolerance;
		}
	}
	
	public double getToleranceThreshold() {
		if (amuTolerance==UNUSED_TOLERANCE) {
			return ppmTolerance;
		} else {
			return amuTolerance;
		}
	}
	
	public double getDeltaScore(double expected, double acquired) {
		if (amuTolerance==UNUSED_TOLERANCE) {
			return (expected-acquired)*1000000.0/expected;
		} else {
			return expected-acquired;
		}
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
		double tolerance=amuTolerance!=UNUSED_TOLERANCE?amuTolerance:Math.max(Math.abs(m1), Math.abs(m2))*percent;
		if (m1+tolerance<m2) return -1;
		if (m1-tolerance>m2) return 1;
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
		if (peaks.length==0) return Optional.empty();
		
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
	public int[] getIndices(TDoubleArrayList peaks, double target) {
		int value=peaks.binarySearch(target);
		// exact match (not likely)
		if (value<0) {
			// insertion point
			value=-(value+1);
		}
		
		TIntArrayList matches=new TIntArrayList();
		// look below
		int index=value;
		while (index>0&&compareTo(peaks.get(index-1), target)==0) {
			matches.add(index-1);
			index--;
		}

		// look up
		index=value;
		while (index<peaks.size()&&compareTo(peaks.get(index), target)==0) {
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
	public int[] getIndices(double[] peaks, double target) {
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
		int[] indicies=getIndices(peaks, target);
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
		
		int targetIndex=0;
		int spectrumIndex=0;
		while (true) {
			double targetMass=targets[targetIndex];
			int compare=compareTo(targetMass, masses[spectrumIndex]);
			if (compare==0) {
				tics[targetIndex]+=intensities[spectrumIndex];
				//libraryIndex++; // could match multiple acquired peaks to the same library peak
				spectrumIndex++;
			} else if (compare>0) {
				spectrumIndex++;
			} else {
				targetIndex++;
			}
			if (targetIndex>=targets.length) break;
			if (spectrumIndex>=masses.length) break;
		}
		
		return tics;
	}
	
	/**
	 * assumes targets and masses are in sorted order (and intensities follows masses)
	 * @param masses
	 * @param intensities
	 * @param targets
	 * @return
	 */
	public Pair<double[], float[]> getIntegratedIntensitiesAndMasses(double[] masses, float[] intensities, double[] targets) {
		float[] tics=new float[targets.length];
		double[] totalMasses=new double[targets.length];
		
		if (targets.length==0||masses.length==0) {
			return new Pair<>(targets, tics);
		}
		
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			double targetMass=targets[libraryIndex];
			int compare=compareTo(targetMass, masses[spectrumIndex]);
			if (compare==0) {
				tics[libraryIndex]+=intensities[spectrumIndex];
				totalMasses[libraryIndex]+=masses[spectrumIndex]*intensities[spectrumIndex];
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
		
		// weighted masses
		for (int i = 0; i < totalMasses.length; i++) {
			if (totalMasses[i]==0.0) {
				totalMasses[i]=targets[i];
			} else {
				totalMasses[i]=totalMasses[i]/tics[i];
			}
		}
		return new Pair<>(totalMasses, tics);
	}
	
	/**
	 * @param peaks -- assumes sorted array of peaks
	 * @param target
	 * @return all matching masses in range
	 */
	public float getIntegratedIntensity(double[] masses, float[] intensities, double target) {
		int[] indicies=getIndices(masses, target);
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
		int[] indicies=getIndices(masses, target);
		float intensity=0.0f;
		for (int i=0; i<indicies.length; i++) {
			if(intensity<intensities[indicies[i]]) {
				intensity=intensities[indicies[i]];
			}
		}
		return intensity;
	}
}
