package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class PecanOneScoringTask extends AbstractPecanScoringTask {
	public PecanOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, TDoubleObjectHashMap<XYPoint> background, PrecursorScanMap precursors,
			int scanAveragingWindow, BlockingQueue<PeptideScoringResult> resultsQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, background, precursors, scanAveragingWindow, resultsQueue, parameters);
	}

	@Override
	protected Nothing process() {
		for (LibraryEntry entry : super.entries) {
			int requiredNumAboveThreshold=(int)(0.5f*entry.getPeptideSeq().length());
			
			int scanAveragingHalfWindow=scanAveragingWindow/2;
			
			float[] rawRTs=new float[super.stripes.size()];
			float[] rawScores=new float[super.stripes.size()];
			float[] bgsubScores=new float[super.stripes.size()];
			float[][] fragmentTraces=new float[entry.getMassArray().length][];
			
			for (int i=0; i<fragmentTraces.length; i++) {
				fragmentTraces[i]=new float[super.stripes.size()];
			}
			
			for (int i=0; i<super.stripes.size(); i++) {
				Stripe stripe=super.stripes.get(i);
				rawRTs[i]=stripe.getScanStartTime();
				
				float[] scores=scorer.getIndividualPeakScores(entry, stripe, false);
				for (int j=0; j<scores.length; j++) {
					fragmentTraces[j][i]=scores[j];
				}
				rawScores[i]+=scorer.score(entry, stripe, precursors);
				
				float rt=stripe.getScanStartTime();
				XYPoint meanStdev=background.get((double)rt);
				if (meanStdev!=null) {
					bgsubScores[i]=(float)(rawScores[i]-meanStdev.x);
				} else {
					bgsubScores[i]=rawScores[i];
				}
			}

			float[] sumRawScores=movingForwardSum(rawScores, scanAveragingWindow);
			float[] sumBgsubScores=movingForwardSum(bgsubScores, scanAveragingWindow);
			float[][] sumFragmentTraces=new float[entry.getIntensityArray().length][];
			for (int i=0; i<sumFragmentTraces.length; i++) {
				sumFragmentTraces[i]=movingForwardSum(fragmentTraces[i], scanAveragingWindow);
			}
			float[] midTime=movingForwardRTAverage(rawRTs, scanAveragingWindow);

			ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
			int[] numAboveThresholdMatches=new int[sumRawScores.length];
			int[] numMatches=new int[sumRawScores.length];
			for (int i=0; i<numAboveThresholdMatches.length; i++) {
				// TODO: this seems questionable that unnormalized intensities are used for individual scores while normalized intensities are used for total scores. -BCS
				float threshold=Math.max(0.0f, sumRawScores[i]/(entry.getPeptideSeq().length()+1));
				for (int j=0; j<sumFragmentTraces.length; j++) {
					if (sumFragmentTraces[j][i]>threshold) {
						numAboveThresholdMatches[i]++;
					}
					if (sumFragmentTraces[j][i]>0.0f) {
						numMatches[i]++;
					}
				}
				
				if (numAboveThresholdMatches[i]>requiredNumAboveThreshold) {
					goodStripes.add(new ScoredIndex(sumBgsubScores[i], i));
				}
			}
			Collections.sort(goodStripes);

			PeptideScoringResult result=new PeptideScoringResult(entry);
			TIntHashSet takenScans=new TIntHashSet();
			int identifiedPeaks=0;
			for (int i=goodStripes.size()-1; i>=0; i--) {
				int index=goodStripes.get(i).y;
				if (takenScans.contains(index)) {
					continue;
				} else {
					int lowerWindow=index-2*scanAveragingWindow; // can't pick anything in twice the peak width
					int upperWindow=index+3*scanAveragingWindow; // +1 to account for the window boundary
					for (int j=lowerWindow; j<upperWindow; j++) {
						takenScans.add(j);
					}
					
					float[][] auxScores=new float[scanAveragingWindow][];
					for (int j=0; j<scanAveragingWindow; j++) {
						Stripe stripe=stripes.get(index+j);
						auxScores[j]=scorer.auxScore(entry, stripe, precursors);
					}
					
					float[] averageAuxScores=new float[auxScores[0].length];
					for (int auxIndex=0; auxIndex<averageAuxScores.length; auxIndex++) {
						for (int scanIndex=0; scanIndex<auxScores.length; scanIndex++) {
							averageAuxScores[auxIndex]+=auxScores[scanIndex][auxIndex];
						}
					}
					for (int j=0; j<averageAuxScores.length; j++) {
						averageAuxScores[j]=averageAuxScores[j]/scanAveragingWindow;
					}

					// averaging forward, so current scan is actually the median for half a window back
					int medianIndex=index+scanAveragingHalfWindow;
					Stripe medianStripe=stripes.get(medianIndex);
					float[] completeAuxArray=General.concatenate(new float[] {numAboveThresholdMatches[index], numMatches[index], midTime[index]}, averageAuxScores, auxScores[scanAveragingHalfWindow]);
					result.addStripe(goodStripes.get(i).x/scanAveragingWindow, completeAuxArray, medianStripe);
					
					if (identifiedPeaks>parameters.getNumberOfReportedPeaks()) {
						// keep N+1 peaks
						break;
					}
					identifiedPeaks++;
				}
			}
			
			TFloatFloatHashMap scoreMap=new TFloatFloatHashMap();
			for (int i=0; i<sumBgsubScores.length; i++) {
				if (numAboveThresholdMatches[i]>=requiredNumAboveThreshold) {
					scoreMap.put(midTime[i], sumBgsubScores[i]);
				} else {
					scoreMap.put(midTime[i], 0.0f);
				}
			}
			
			result.setTrace(new XYTrace(scoreMap, GraphType.line, entry.getPeptideModSeq()));
			
			try {
				resultsQueue.put(result);
			} catch (InterruptedException ie) {
				Logger.errorLine("Analysis interrupted!");
				Logger.errorException(ie);
				return Nothing.NOTHING;
			}
		}
		return Nothing.NOTHING;
	}
}
