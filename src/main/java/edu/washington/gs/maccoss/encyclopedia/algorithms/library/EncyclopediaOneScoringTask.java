package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.set.hash.TIntHashSet;

public class EncyclopediaOneScoringTask extends AbstractLibraryScoringTask {
	
	public EncyclopediaOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue,
			SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
	}
	
	private static final int movingAverageLength=10;
	private static final int peaksKept=5;

	@Override
	protected Nothing process() {
		for (LibraryEntry entry : super.entries) {
			PeptideScoringResult result=new PeptideScoringResult(entry);
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
			
			float[][] scores=new float[super.stripes.size()][];
			float[] primary=new float[super.stripes.size()];
			for (int i=0; i<super.stripes.size(); i++) {
				Stripe stripe=super.stripes.get(i);
				scores[i]=scorer.auxScore(entry, stripe, predictedIsotopeDistribution, precursors);
				primary[i]=scores[i][0];
				//primary[i]=scorer.score(entry, stripe, predictedIsotopeDistribution, precursors);
			}
			
			float[] averagePrimary=movingCenteredAverage(primary, movingAverageLength);

			ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
			for (int i=0; i<averagePrimary.length; i++) {
				goodStripes.add(new ScoredIndex(primary[i], i));
			}
			Collections.sort(goodStripes);

			TIntHashSet takenScans=new TIntHashSet();
			int identifiedPeaks=0;
			for (int i=goodStripes.size()-1; i>=0; i--) {
				float score=goodStripes.get(i).x;
				int index=goodStripes.get(i).y;
				if (takenScans.contains(index)) {
					continue;
					
				} else {
					float[] auxScoreArray=scores[index];
					result.addStripe(score, auxScoreArray, super.stripes.get(index));
					
					// block out a 40 scan window
					int lowerWindow=index-2*movingAverageLength;
					int upperWindow=index+2*movingAverageLength;
					for (int j=lowerWindow; j<upperWindow; j++) {
						takenScans.add(j);
					}
					
					if (identifiedPeaks>peaksKept) {
						// keep N+1 peaks
						break;
					}
					identifiedPeaks++;
				}
			}
			
			resultsQueue.add(result);
		}
		return Nothing.NOTHING;
	}

	
}
