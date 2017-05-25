package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.list.array.TDoubleArrayList;

public class CorrelationPeakFrequencyCalculator {
	private final MassTolerance tolerance;
	private final ArrayList<Count> peaks=new ArrayList<Count>();
	
	public CorrelationPeakFrequencyCalculator(MassTolerance tolerance) {
		this.tolerance=tolerance;
	}
	
	public ArrayList<XYPoint> toPoints() {
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		for (Count peak : peaks) {
			points.add(new XYPoint(peak.mass, peak.count));
		}
		return points;
	}

	public double[] getTopNMasses(int n) {
		ArrayList<Count> clone=new ArrayList<CorrelationPeakFrequencyCalculator.Count>(peaks);
		Collections.sort(clone);
		
		for (Count count : clone) {
			System.out.println(count.count+"\t"+count.mass+"\t"+count.maxIntensity);
		}
		
		TDoubleArrayList masses=new TDoubleArrayList();
		int count=0;
		for (int i=clone.size()-1; i>=0; i--) {
			masses.add(clone.get(i).mass);
			count++;
			if (count>=n) break;
		}
		return masses.toArray();
	}

	public void increment(double target, float intensity, float correlation) {
		int value=binarySearch(peaks, target);
		if (value>=0) {
			peaks.get(value).increment(intensity, correlation);
		} else {
			value=-(value+1);
			peaks.add(value, new Count(target, intensity, correlation));
		}
	}

	/*
	 * NOTE: Does not compare with Comparable<Count> (based on the count), but instead uses MassTolerance (based on the mass)
	 * Design taken from:
	 * Arrays.binarySearch0(double[] a, int fromIndex, int toIndex, double key)
	 */
	private int binarySearch(ArrayList<Count> peaks, double key) {
		int low=0;
		int high=peaks.size()-1;

		while (low<=high) {
			int mid=(low+high)>>>1;
			double midVal=peaks.get(mid).mass;
			
			if (tolerance.compareTo(midVal, key)==0) {
				return mid; // within tolerance match
			} else if (midVal<key) {
				low=mid+1; // Neither val is NaN, thisVal is smaller
			} else if (midVal>key) {
				high=mid-1; // Neither val is NaN, thisVal is larger
			} else {
				throw new EncyclopediaException("Mass tolerance did not match identical double match (or there are NaNs)!");
			}
		}
		return -(low+1); // key not found.
	}
	
	public class Count implements Comparable<Count> {
		private final double mass;
		private float maxIntensity=0.0f;
		private float count=0.0f;
		public Count(double mass, float intensity, float correlation) {
			this.mass=mass;
			this.maxIntensity=intensity;
			count=correlation;
		}
		public float increment(float intensity, float correlation) {
			if (intensity>maxIntensity) maxIntensity=intensity;
			count=count+correlation;
			return count;
		}
		
		@Override
		/**
		 * sorts on count first then on mass so that the highest count, highest intensity peaks are at the end
		 */
		public int compareTo(Count o) {
			if (o==null) return 1;
			int c=Float.compare(count, o.count);
			if (c!=0) return c;
			return Double.compare(maxIntensity, o.maxIntensity);
		}
	}
}
