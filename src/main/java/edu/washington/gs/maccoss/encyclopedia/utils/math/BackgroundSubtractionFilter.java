package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;

public class BackgroundSubtractionFilter {

	
	public static float[] movingForwardRTAverage(float[] rts, int scanAveragingWindow) {
		// like moving sum, this approach drops the first scanAveragingWindow-1 scans
		float[] avgRTs=new float[rts.length-scanAveragingWindow];
		for (int i=0; i<avgRTs.length; i++) {
			avgRTs[i]=(rts[i]+rts[i+scanAveragingWindow-1])/2.0f;
		}
		return avgRTs;
	}
	
	public static float[] movingForwardSum(float[] scores, int scanAveragingWindow) {
		// moving sum, this approach drops the first scanAveragingWindow-1 scans
		float[] sumScores=new float[scores.length-scanAveragingWindow];
		for (int i=0; i<sumScores.length; i++) {
			for (int j=0; j<scanAveragingWindow; j++) {
				sumScores[i]+=scores[i+j];
			}
		}
		return sumScores;
	}
	
	public static float[] movingForwardAverage(float[] scores, int scanAveragingWindow) {
		// moving sum, this approach drops the first scanAveragingWindow-1 scans
		float[] sumScores=new float[scores.length-scanAveragingWindow];
		for (int i=0; i<sumScores.length; i++) {
			for (int j=0; j<scanAveragingWindow; j++) {
				sumScores[i]+=scores[i+j];
			}
			sumScores[i]=sumScores[i]/scanAveragingWindow;
		}
		return sumScores;
	}

	public static float[] fastMovingCenteredAverage(float[] scores, int scanAveragingWindow) {
		// bad params cases
		if (scanAveragingWindow<=1) {
			return scores;
		} else if (scores.length<=scanAveragingWindow) {
			float average=General.mean(scores);
			float[] newScores=new float[scores.length];
			Arrays.fill(newScores, average);
			return newScores;
		}
		
		// generate the average of the first window 
		float windowSum=0.0f;
		for (int i = 0; i < scanAveragingWindow; i++) {
			windowSum+=scores[i];
		}
		
		// populate before the window with the left edge
		float[] sumScores=new float[scores.length];
		int scanAveragingMargin;
		if (General.isEven(scanAveragingWindow)) {
			scanAveragingMargin=scanAveragingWindow/2;
		} else {
			scanAveragingMargin=(scanAveragingWindow+1)/2;
		}
		for (int i = 0; i <= scanAveragingMargin; i++) {
			sumScores[i]=windowSum;
		}
		
		// then add from the new value / subtract from the first value in the window
		for (int i = scanAveragingWindow+1; i < scores.length; i++) {
			float first=scores[i-1-scanAveragingWindow];
			float newValueToAdd=scores[i];
			windowSum=windowSum-first+newValueToAdd;
			sumScores[i-scanAveragingMargin]=windowSum;
		}

		// populate after the window with the right edge
		for (int i = scores.length-scanAveragingMargin; i < scores.length; i++) {
			sumScores[i]=windowSum;
		}
		
		// return the average
		return General.divide(sumScores, (float)scanAveragingWindow);
	}
	
	/**
	 * algorithm maintains a sorted window loosly following ideas in Mohanty 2003 (https://dcc-backup.ligo.org/public/0027/T030168/000/T030168-00.pdf)
	 */
	public static float[] fastMovingMedian(float[] scores, int scanAveragingWindow) {
		// bad params cases
		if (scanAveragingWindow<=1) {
			return scores;
		} else if (scores.length<=scanAveragingWindow) {
			float average=General.mean(scores);
			float[] newScores=new float[scores.length];
			Arrays.fill(newScores, average);
			return newScores;
		}
		
		float[] window=new float[scanAveragingWindow];
		System.arraycopy(scores, 0, window, 0, scanAveragingWindow);
		Arrays.sort(window);
		float median=getMedianOfSortedArray(window);

		// populate before the window with the left edge
		float[] medianScores=new float[scores.length];
		int scanAveragingMargin;
		if (General.isEven(scanAveragingWindow)) {
			scanAveragingMargin=scanAveragingWindow/2;
		} else {
			scanAveragingMargin=(scanAveragingWindow+1)/2;
		}
		for (int i = 0; i <= scanAveragingMargin; i++) {
			medianScores[i]=median;
		}
		
		// then add from the new value / subtract from the first value in the window
		int windowLengthMinusOne = window.length-1;
		for (int i = scanAveragingWindow; i < scores.length; i++) {
			float first=scores[i-scanAveragingWindow];
			float newValueToAdd=scores[i];
			
			// find the location of the old value and swap with the new value
			int index=Arrays.binarySearch(window, first);
			if (index<0) throw new EncyclopediaException("Median calculation error with missing value "+first+"("+i+","+(i-scanAveragingWindow)+"), index="+index+", ");
			window[index]=newValueToAdd;
			
			// twiddle up or down to maintain the sort order
			if (newValueToAdd>first) {
				// move up
				while (index<windowLengthMinusOne&&window[index]>window[index+1]) {
					float temp=window[index+1];
					window[index+1]=window[index];
					window[index]=temp;
					index=index+1;
				}
			} else if (newValueToAdd<first) {
				// move down
				while (index>0&&window[index]<window[index-1]) {
					float temp=window[index-1];
					window[index-1]=window[index];
					window[index]=temp;
					index=index-1;
				}
			}

			median=getMedianOfSortedArray(window);
			medianScores[i-scanAveragingMargin+1]=median;
		}

		// populate after the window with the right edge
		for (int i = scores.length-scanAveragingMargin; i < scores.length; i++) {
			medianScores[i]=median;
		}
		
		return medianScores;
	}
	
	private static float getMedianOfSortedArray(float[] scores) {
		if (General.isEven(scores.length)) {
			// if 4, we want index 1&2 (the 2nd and 3rd values) 
			return (scores[scores.length/2-1]+scores[scores.length/2])/2.0f;
		} else {
			// if 3, we want index 1 (the 2nd value)
			return scores[scores.length/2];
		}
	}
	
	public static float[] backgroundSubtractMovingAverage(float[] scores, int scanAveragingWindow) {
		float[] movingAverage=fastMovingCenteredAverage(scores, scanAveragingWindow);
		float[] subtract=General.subtract(scores, movingAverage);
		for (int i = 0; i < subtract.length; i++) {
			if (subtract[i]<0.0f) subtract[i]=0.0f;
		}
		return subtract;
	}
	/**
	 * assumes even incremented X values!
	 * 
	 * @param trace
	 * @return
	 */
	public static XYTrace backgroundSubtractMovingMedian(XYTraceInterface trace, int scanAveragingWindow) {
		Pair<double[], double[]> values=trace.toArrays();
		float[] smoothedY=backgroundSubtractMovingMedian(General.toFloatArray(values.y), scanAveragingWindow);
		return new XYTrace(values.x, General.toDoubleArray(smoothedY), trace.getType(), trace.getName(), trace.getColor().orElse(null), trace.getThickness().orElse(null));
	}
	
	public static float[] backgroundSubtractMovingMedian(float[] scores, int scanAveragingWindow) {
		float[] movingAverage=fastMovingMedian(scores, scanAveragingWindow);
		float[] subtract=General.subtract(scores, movingAverage);
		for (int i = 0; i < subtract.length; i++) {
			if (subtract[i]<0.0f) subtract[i]=0.0f;
		}
		return subtract;
	}
	
	public static float[] movingCenteredAverage(float[] scores, int scanAveragingWindow) {
		// moving sum on background subtracted scores, this approach uses less data for the first and last scanAveragingMargin scans
		int scanAveragingMargin=(scanAveragingWindow-1)/2;
		
		float[] sumScores=new float[scores.length];
		for (int i=0; i<scores.length; i++) {
			float sum=0.0f;
			for (int j=0; j<scanAveragingWindow; j++) {
				int index=i+j-scanAveragingMargin;
				if (index>=0&&index<scores.length) {
					sum+=scores[index];
				}
			}
			sumScores[i]=sum/scanAveragingWindow;
		}
		return sumScores;
	}
	
	public static float[] gaussianCenteredAverage(float[] scores, int scanAveragingWindow) {
		float mean=(scanAveragingWindow-1)/2.0f; // -1 to get the real center index
		float stdev=(scanAveragingWindow-1)/6.0f; // minus the center, calculate the stdev (6 assumes a peak is 3 stdevs on either side of the center)
		Gaussian g=new Gaussian(mean, stdev, 1.0f);
		float[] probs=new float[scanAveragingWindow];
		for (int i=0; i<scanAveragingWindow; i++) {
			probs[i]=(float)g.getPDF(i);
		}
		
		// moving sum on background subtracted scores, this approach uses less data for the first and last scanAveragingMargin scans
		int scanAveragingMargin=(scanAveragingWindow-1)/2;
		
		float[] sumScores=new float[scores.length];
		for (int i=0; i<scores.length; i++) {
			float sum=0.0f;
			for (int j=0; j<scanAveragingWindow; j++) {
				float prob=probs[j];
				int index=i+j-scanAveragingMargin;
				if (index>=0&&index<scores.length) {
					sum+=prob*scores[index];
				}
			}
			sumScores[i]=sum;
		}
		return sumScores;
	}

}
