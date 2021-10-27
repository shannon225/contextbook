package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
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
	
	public static float[] movingCenteredSum(float[] scores, int scanAveragingWindow) {
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
			sumScores[i]=sum;
		}
		return sumScores;
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
