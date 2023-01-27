package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Gaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;

public abstract class AbstractLibraryScoringTask extends ThreadableTask<Nothing> {
	/**
	 * must be immutable!
	 */
	protected final PSMScorer scorer;
	protected final ArrayList<LibraryEntry> entries;
	protected final ArrayList<FragmentScan> stripes;
	protected final PrecursorScanMap precursors;
	protected final BlockingQueue<AbstractScoringResult> resultsQueue;
	protected final SearchParameters parameters;

	/**
	 * scorer must be a 
	 * @param scorer
	 * @param entries
	 * @param stripes
	 * @param background
	 * @param precursors
	 * @param scanAveragingMargin
	 */
	public AbstractLibraryScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue, SearchParameters parameters) {
		this.scorer=scorer;
		this.entries=entries;
		this.stripes=stripes;
		this.precursors=precursors;
		this.resultsQueue=resultsQueue;
		this.parameters=parameters;
	}

	@Override
	public String getTaskName() {
		StringBuilder sb=new StringBuilder();
		for (LibraryEntry entry : entries) {
			if (sb.length()>0) {
				sb.append(',');
			}
			sb.append(entry.getPeptideModSeq());
		}
		return sb.toString();
	}
	
	protected float[] movingForwardRTAverage(float[] rts, int scanAveragingWindow) {
		// like moving sum, this approach drops the first scanAveragingWindow-1 scans
		float[] avgRTs=new float[rts.length-scanAveragingWindow];
		for (int i=0; i<avgRTs.length; i++) {
			avgRTs[i]=(rts[i]+rts[i+scanAveragingWindow-1])/2.0f;
		}
		return avgRTs;
	}
	
	protected float[] movingForwardSum(float[] scores, int scanAveragingWindow) {
		// moving sum, this approach drops the first scanAveragingWindow-1 scans
		float[] sumScores=new float[scores.length-scanAveragingWindow];
		for (int i=0; i<sumScores.length; i++) {
			for (int j=0; j<scanAveragingWindow; j++) {
				sumScores[i]+=scores[i+j];
			}
		}
		return sumScores;
	}
	
	protected float[] movingForwardAverage(float[] scores, int scanAveragingWindow) {
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
	public static float[] fastMovingQuickMedian(float[] scores, int scanAveragingWindow) {
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
		for (int i = scanAveragingWindow; i < scores.length; i++) {
			float first=scores[i-scanAveragingWindow];
			float newValueToAdd=scores[i];
			
			//System.out.println((i-scanAveragingWindow)+"\tSwap "+first+" for "+newValueToAdd+" --> "+General.toString(window));
			
			int index=Arrays.binarySearch(window, first);
			if (index<0) throw new EncyclopediaException("Median calculation error with missing value "+first+"("+i+","+(i-scanAveragingWindow)+"), index="+index+", ");
			window[index]=newValueToAdd;
			Arrays.sort(window);

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
	
	public static float[] backgroundSubtractMovingMedian(float[] scores, int scanAveragingWindow) {
		float[] movingAverage=fastMovingQuickMedian(scores, scanAveragingWindow);
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
