package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.IndexedObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class PecanScoringTask extends PeptideScoringTask {
	private final int scanAveragingMargin;
	protected final TDoubleObjectHashMap<XYPoint> background; // if not null, then score using zscore (otherwise use raw score)

	public PecanScoringTask(PecanRawScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, TDoubleObjectHashMap<XYPoint> background, int scanAveragingMargin) {
		super(scorer, entries, stripes);
		this.background=background;
		this.scanAveragingMargin=scanAveragingMargin;
	}
	private PecanRawScorer getScorer() {
		return (PecanRawScorer)super.scorer;
	}

	@Override
	protected HashMap<LibraryEntry, PeptideScoringResult> process() {
		HashMap<LibraryEntry, PeptideScoringResult> map=new HashMap<LibraryEntry, PeptideScoringResult>();
		
		for (LibraryEntry entry : super.entries) {
			float requiredNumAboveThreshold=0.5f*entry.getPeptideSeq().length();
			
			TFloatArrayList rawScores=new TFloatArrayList();
			TFloatArrayList backgroundSubtractedScores=new TFloatArrayList();
			TFloatArrayList numAboveThresholdMatches=new TFloatArrayList();
			
			for (Stripe stripe : super.stripes) {
				float[] scores=getScorer().subScore(entry, stripe);
				float rawScore=scores[0];
				float numAboveThreshold=scores[1];
				
				rawScores.add(rawScore);
				numAboveThresholdMatches.add(numAboveThreshold);

				float rt=stripe.getScanStartTime();
				XYPoint meanStdev=background.get((double)rt);
				if (meanStdev!=null) {
					backgroundSubtractedScores.add((float)(rawScore-meanStdev.x));
				} else {
					backgroundSubtractedScores.add(rawScore);
				}
			}

			ArrayList<ScoredObject<IndexedObject<Stripe>>> goodStripes=new ArrayList<ScoredObject<IndexedObject<Stripe>>>();

			TFloatArrayList windowedBackgroundSubtractedScores=new TFloatArrayList();
			int scanAveragingWindow=2*scanAveragingMargin+1;
			int scanExcludingWindow=2*scanAveragingWindow+1;
			// moving average on background subtracted scores, this approach uses less data for the first and last scanAveragingMargin scans
			TIntFloatHashMap scoreByIndex=new TIntFloatHashMap();
			for (int i=0; i<backgroundSubtractedScores.size(); i++) {
				float sum=0.0f;
				int count=0;
				for (int j=0; j<scanAveragingWindow; j++) {
					int index=i+j-scanAveragingMargin;
					if (index>=0&&index<backgroundSubtractedScores.size()) {
						sum+=backgroundSubtractedScores.get(index);
						count++;
					}
				}
				float smoothedScore=sum/count;
				windowedBackgroundSubtractedScores.add(smoothedScore);
				scoreByIndex.put(i, smoothedScore);

				if (numAboveThresholdMatches.get(i)>=requiredNumAboveThreshold) {
					goodStripes.add(new ScoredObject<IndexedObject<Stripe>>(smoothedScore, new IndexedObject<Stripe>(i, stripes.get(i))));
				}
			}
			Collections.sort(goodStripes);

			PeptideScoringResult result=new PeptideScoringResult();
			TIntHashSet takenScans=new TIntHashSet();
			for (int i=goodStripes.size()-1; i>=0; i--) {
				IndexedObject<Stripe> stripe=goodStripes.get(i).y;
				if (takenScans.contains(stripe.x)) {
					continue;
				} else {
					float[] averageAuxScores=null;
					float total=0.0f;
					for (int j=0; j<scanExcludingWindow; j++) {
						int index=stripe.x-scanAveragingWindow+j;
						float indexScore=scoreByIndex.get(index);
						
						takenScans.add(index);
						float[] auxScores=getScorer().auxScore(entry, stripe.y);
						if (indexScore>0) {
							if (averageAuxScores==null) {
								averageAuxScores=General.multiply(auxScores, indexScore);
							} else {
								averageAuxScores=General.add(averageAuxScores, General.multiply(auxScores, indexScore));
							}
						}
						total+=indexScore;
					}
					if (averageAuxScores!=null&&total>=0.0f) {
						averageAuxScores=General.multiply(averageAuxScores, 1.0f/total);
					} else {
						averageAuxScores=getScorer().getAuxScorer().getMissingDataScores();
					}
					
					result.addStripe(goodStripes.get(i).x, averageAuxScores, stripe.y);
				}
			}
			
			//EValueCalculator calculator=new EValueCalculator(scoreMap);
			TFloatFloatHashMap scoreMap=new TFloatFloatHashMap();
			for (int i=0; i<super.stripes.size(); i++) {
				if (numAboveThresholdMatches.get(i)>=requiredNumAboveThreshold) {
					scoreMap.put(super.stripes.get(i).getScanStartTime(), windowedBackgroundSubtractedScores.get(i));
				} else {
					scoreMap.put(super.stripes.get(i).getScanStartTime(), 0.0f);
				}
			}
			
			result.setTrace(new XYTrace(scoreMap, GraphType.line, entry.getPeptideModSeq()));
			map.put(entry, result);
		}
		return map;
	}
}
