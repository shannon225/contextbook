package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class Correlation {
	public static float dotProduct(Spectrum predicted, Spectrum acquired, MassTolerance tolerance) {
		double[] libraryMasses=predicted.getMassArray();
		float[] libraryIntensities=predicted.getIntensityArray();
		
		double[] spectrumMasses=acquired.getMassArray();
		float[] spectrumIntensities=acquired.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return 0.0f;
		
		float sum=0.0f;
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			int compare=tolerance.compareTo(libraryMasses[libraryIndex], spectrumMasses[spectrumIndex]);
			if (compare==0) {
				sum+=libraryIntensities[libraryIndex]*spectrumIntensities[spectrumIndex];
				libraryIndex++;
				spectrumIndex++;
			} else if (compare>0) {
				spectrumIndex++;
			} else {
				libraryIndex++;
			}
			if (libraryIndex>=libraryMasses.length) break;
			if (spectrumIndex>=spectrumMasses.length) break;
		}
		
		return sum;
	}
	
	public static double getSpearmans(Spectrum predicted, Spectrum acquired, MassTolerance tolerance) {
		Pair<double[], double[]> arrays=getArrays(predicted, acquired, tolerance);
		if (arrays==null) return 0.0;
		return getSpearmans(arrays.x, arrays.y);
	}
	
	public static double getPearsons(Spectrum predicted, Spectrum acquired, MassTolerance tolerance) {
		Pair<double[], double[]> arrays=getArrays(predicted, acquired, tolerance);
		if (arrays==null) return 0.0;
		return getPearsons(arrays.x, arrays.y);
	}

	private static Pair<double[], double[]> getArrays(Spectrum predicted, Spectrum acquired, MassTolerance tolerance) {
		double[] libraryMasses=predicted.getMassArray();
		float[] libraryIntensities=predicted.getIntensityArray();
		
		double[] spectrumMasses=acquired.getMassArray();
		float[] spectrumIntensities=acquired.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return null;
		
		int libraryIndex=0;
		int spectrumIndex=0;
		ArrayList<XYPoint> points=new ArrayList<>();
		while (true) {
			int compare=tolerance.compareTo(libraryMasses[libraryIndex], spectrumMasses[spectrumIndex]);
			if (compare==0) {
				points.add(new XYPoint(libraryIntensities[libraryIndex], spectrumIntensities[spectrumIndex]));
				libraryIndex++;
				spectrumIndex++;
			} else if (compare>0) {
				points.add(new XYPoint(0.0, spectrumIntensities[spectrumIndex]));
				spectrumIndex++;
			} else {
				points.add(new XYPoint(libraryIntensities[libraryIndex], 0.0));
				libraryIndex++;
			}
			if (libraryIndex>=libraryMasses.length) {
				for (int i=spectrumIndex; i<spectrumIntensities.length; i++) {
					points.add(new XYPoint(0.0, spectrumIntensities[spectrumIndex]));
				}
				break;
			}
			if (spectrumIndex>=spectrumMasses.length) {

				for (int i=libraryIndex; i<libraryIntensities.length; i++) {
					points.add(new XYPoint(0.0, libraryIntensities[libraryIndex]));
				}
				break;
			}
		}
		
		Pair<double[], double[]> arrays=XYTrace.toArrays(points);
		return arrays;
	}
	
	public static double getSpearmans(double[] x, double[] y) {
		return getPearsons(rank(x), rank(y));
	}
	
	public static double getPearsons(double[] x, double[] y) {
		
		double xBar=General.mean(x);
		double yBar=General.mean(y);
		
		double numerator=0.0;
		double xSS=0.0;
		double ySS=0.0;
		for (int i=0; i<y.length; i++) {
			double xDiff=x[i]-xBar;
			double yDiff=y[i]-yBar;
			numerator+=xDiff*yDiff;
			xSS+=xDiff*xDiff;
			ySS+=yDiff*yDiff;
		}
		if (xSS==0||ySS==0) {
			return 0.0;
		}
		return numerator/Math.sqrt(xSS*ySS);
	}
	
	private static double[] rank(double[] values) {
		double[] sorted=values.clone();
		Arrays.sort(sorted);
		double[] sortedRanks=new double[sorted.length];

		// basic rank
		for (int i=0; i<sorted.length; i++) {
			sortedRanks[i]=i+1;
		}

		for (int i=0; i<sorted.length; i++) {
			int start=i;
			int stop=i;

			// find ties
			boolean ties=false;
			while (++stop<sorted.length&&sorted[start]==sorted[stop]) {
				ties=true;
			}

			// substitute rank average for ties
			if (stop-start>1&&ties) {
				double avg=0;
				for (int j=start; j<stop; j++) {
					avg+=sortedRanks[j];
				}
				avg=avg/(stop-start);

				for (int x=start; x<stop; x++) {
					sortedRanks[x]=avg;
				}
			}
			
			// advance i to end of stop
			i=stop-1;
		}
		double[] ranksInOrder=new double[sortedRanks.length];
		for (int i=0; i<sortedRanks.length; i++) {
			int index=Arrays.binarySearch(sorted, values[i]);
			ranksInOrder[i]=sortedRanks[index];
		}

		return ranksInOrder;
	}
}
