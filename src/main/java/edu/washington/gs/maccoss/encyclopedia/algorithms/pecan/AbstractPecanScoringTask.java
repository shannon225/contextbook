package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringTask;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public class AbstractPecanScoringTask extends PeptideScoringTask {
	protected final int scanAveragingMargin;
	protected final TDoubleObjectHashMap<XYPoint> background; // if not null, then score using zscore (otherwise use raw score)

	/**
	 * scorer must be a 
	 * @param scorer
	 * @param entries
	 * @param stripes
	 * @param background
	 * @param precursors
	 * @param scanAveragingMargin
	 */
	public AbstractPecanScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, TDoubleObjectHashMap<XYPoint> background, PrecursorScanMap precursors, int scanAveragingMargin) {
		super(scorer, entries, stripes, precursors);
		this.background=background;
		this.scanAveragingMargin=scanAveragingMargin;
	}
	
	protected float[] movingSum(float[] scores, int scanAveragingWindow) {
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
