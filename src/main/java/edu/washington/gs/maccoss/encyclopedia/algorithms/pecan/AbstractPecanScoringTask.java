package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public abstract class AbstractPecanScoringTask extends ThreadableTask<Nothing> {
	/**
	 * must be immutable!
	 */
	protected final PSMScorer scorer;
	protected final ArrayList<LibraryEntry> entries;
	protected final ArrayList<Stripe> stripes;
	protected final PrecursorScanMap precursors;
	protected final int scanAveragingWindow;
	protected final TDoubleObjectHashMap<XYPoint> background; // if not null, then score using zscore (otherwise use raw score)
	protected final BlockingQueue<PeptideScoringResult> resultsQueue;
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
	public AbstractPecanScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, TDoubleObjectHashMap<XYPoint> background, PrecursorScanMap precursors, int scanAveragingWindow, BlockingQueue<PeptideScoringResult> resultsQueue, SearchParameters parameters) {
		this.scorer=scorer;
		this.entries=entries;
		this.stripes=stripes;
		this.precursors=precursors;
		this.background=background;
		this.scanAveragingWindow=scanAveragingWindow;
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
	
	protected float[] movingCenteredSum(float[] scores, int scanAveragingWindow) {
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
}
