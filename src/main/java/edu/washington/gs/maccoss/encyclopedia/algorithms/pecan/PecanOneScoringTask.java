package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.SortLaterXYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BackgroundSubtractionFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class PecanOneScoringTask extends AbstractPecanScoringTask {

	public PecanOneScoringTask(PSMPeakScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, TDoubleObjectHashMap<XYPoint>[] background, PrecursorScanMap precursors,
			int scanAveragingWindow, BlockingQueue<AbstractScoringResult> resultsQueue, PecanSearchParameters parameters) {
		super(scorer, entries, stripes, background, precursors, scanAveragingWindow, resultsQueue, parameters);
	}

	@Override
	protected Nothing process() {
		for (LibraryEntry entry : super.entries) {
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());

			int requiredNumAboveThreshold=(int)(getPecanSearchParameters().getBeta()*entry.getMassArray().length);

			int scanAveragingHalfWindow=scanAveragingWindow/2;

			float[] rawRTs=new float[super.stripes.size()];
			float[] rawScores=new float[super.stripes.size()];
			float[] bgsubScores=new float[super.stripes.size()];
			float[] backgroundScores=new float[super.stripes.size()];
			float[] zScores=new float[super.stripes.size()];
			float[][] fragmentTraces=new float[entry.getMassArray().length][];
			float[][] fragmentDeltaMasses=new float[entry.getMassArray().length][];

			for (int i=0; i<fragmentTraces.length; i++) {
				fragmentTraces[i]=new float[super.stripes.size()];
				fragmentDeltaMasses[i]=new float[super.stripes.size()];
			}
			
			TDoubleObjectHashMap<XYPoint> chargeStateBackground;
			if (background!=null&&background.length>entry.getPrecursorCharge()) {
				chargeStateBackground=background[entry.getPrecursorCharge()-1];
			} else {
				chargeStateBackground=null;
			}

			for (int i=0; i<super.stripes.size(); i++) {
				FragmentScan stripe=super.stripes.get(i);
				rawRTs[i]=stripe.getScanStartTime();

				// known to be PSMPeakScorer because of constructor
				PeakScores[] scores=((PSMPeakScorer)scorer).getIndividualPeakScores(entry, stripe, true);
				for (int j=0; j<scores.length; j++) {
					if (scores[j]!=null) {
						fragmentTraces[j][i]=scores[j].getScore();
						fragmentDeltaMasses[j][i]=scores[j].getDeltaMass();
					} else {
						fragmentTraces[j][i]=0.0f;
						fragmentDeltaMasses[j][i]=0.0f;
					}
				}
				rawScores[i]+=PeakScores.sumScores(scores); // dot product

				float rt=stripe.getScanStartTime();
				if (chargeStateBackground!=null) {
					XYPoint meanStdev=chargeStateBackground.get((double)rt);
					if (meanStdev!=null) {
						backgroundScores[i]=(float)meanStdev.x;
						bgsubScores[i]=(float)(rawScores[i]-meanStdev.x);
						if (meanStdev.y==0) {
							zScores[i]=0.0f;
						} else {
							zScores[i]=(float)((rawScores[i]-meanStdev.x)/meanStdev.y);
						}
					} else {
						bgsubScores[i]=rawScores[i];
						zScores[i]=0.0f;
					}
				} else {
					bgsubScores[i]=rawScores[i];
					zScores[i]=0.0f;
				}
			}

			float[] sumRawScores=BackgroundSubtractionFilter.movingForwardSum(rawScores, scanAveragingWindow);
			float[] sumBgsubScores=BackgroundSubtractionFilter.movingForwardSum(bgsubScores, scanAveragingWindow);
			float[] sumZScores=BackgroundSubtractionFilter.movingForwardSum(zScores, scanAveragingWindow);
			float[][] sumFragmentTraces=new float[entry.getIntensityArray().length][];
			float[][] sumFragmentDeltaMasses=new float[entry.getIntensityArray().length][];
			for (int i=0; i<sumFragmentTraces.length; i++) {
				sumFragmentTraces[i]=BackgroundSubtractionFilter.movingForwardSum(fragmentTraces[i], scanAveragingWindow);
				sumFragmentDeltaMasses[i]=BackgroundSubtractionFilter.movingForwardAverage(fragmentDeltaMasses[i], scanAveragingWindow);
			}
			float[] midTime=BackgroundSubtractionFilter.movingForwardRTAverage(rawRTs, scanAveragingWindow);

			ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
			int[] numAboveThresholdMatches=new int[sumRawScores.length];
			int[] numMatches=new int[sumRawScores.length]; // calculated but not reported
			
			float scoreThresholdHyperParameter=(float)Math.pow(entry.getMassArray().length, getPecanSearchParameters().getAlpha());
			for (int i=0; i<numAboveThresholdMatches.length; i++) {
				float threshold=Math.max(Float.MIN_VALUE, sumRawScores[i]/scoreThresholdHyperParameter);
				for (int j=0; j<sumFragmentTraces.length; j++) {
					if (sumFragmentTraces[j][i]>=threshold) {
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

			float[] sortedBGScores=sumBgsubScores.clone();
			Arrays.sort(sortedBGScores);

			AbstractScoringResult result=new PeptideScoringResult(entry);
			TIntHashSet takenScans=new TIntHashSet();
			int identifiedPeaks=0;
			for (int i=goodStripes.size()-1; i>=0; i--) {
				int index=goodStripes.get(i).y;
				if (takenScans.contains(index)) {
					continue;

				} else {
					float endScanTime=stripes.get(index+scanAveragingWindow-1).getScanStartTime();
					float startScanTime=stripes.get(index).getScanStartTime();
					float duration=endScanTime-startScanTime;

					// one delta mass per fragment ion
					float[] weightedFragmentDeltaMass=new float[sumFragmentDeltaMasses.length];
					for (int d=0; d<sumFragmentDeltaMasses.length; d++) {
						float totalSumRawScore=0.0f;
						for (int j=0; j<scanAveragingWindow; j++) {
							int scanIndex=index+j;
							if (scanIndex<sumRawScores.length) {
								totalSumRawScore+=sumRawScores[scanIndex];
								weightedFragmentDeltaMass[d]+=sumRawScores[scanIndex]*fragmentDeltaMasses[d][scanIndex];
							}
						}
						if (totalSumRawScore>0) {
							weightedFragmentDeltaMass[d]=weightedFragmentDeltaMass[d]/totalSumRawScore;
						}
					}
					float fragmentDeltaMassAverage=General.mean(weightedFragmentDeltaMass);
					float fragmentDeltaMassVariance=General.variance(weightedFragmentDeltaMass);

					float[][] auxScores=new float[scanAveragingWindow][];
					for (int j=0; j<scanAveragingWindow; j++) {
						FragmentScan stripe=stripes.get(index+j);
						auxScores[j]=scorer.auxScore(entry, stripe, predictedIsotopeDistribution, precursors);
					}

					float[] averageAuxScores=new float[auxScores[0].length];
					for (int auxIndex=0; auxIndex<averageAuxScores.length; auxIndex++) {
						for (int scanIndex=0; scanIndex<auxScores.length; scanIndex++) {
							if (auxScores[scanIndex].length>auxIndex) {
								averageAuxScores[auxIndex]+=auxScores[scanIndex][auxIndex];
							}
						}
					}

					for (int j=0; j<averageAuxScores.length; j++) {
						averageAuxScores[j]=averageAuxScores[j]/scanAveragingWindow;
					}

					// TODO indexing these is hokey, and should be more firmly rooted in the scoring system
					float maxIDP=0.0f; // IDP is the last score of the PecanAuxillaryScorer
					int midIndex=auxScores.length/2;
					float midIDP=auxScores[midIndex][auxScores[midIndex].length-1];
					float precursorPPMVariance=0.0f; // PPM is the second to
														// last score of
														// PecanAuxillaryScorer
					for (int scanIndex=0; scanIndex<auxScores.length; scanIndex++) {
						// TODO indexing these is hokey, and should be more firmly rooted in the scoring system
						int idpIndex=auxScores[scanIndex].length-1;
						int ppmIndex=auxScores[scanIndex].length-2;

						if (auxScores[scanIndex][idpIndex]>maxIDP) {
							maxIDP=auxScores[scanIndex][idpIndex];
						}

						// FIXME VARIANCE NEEDS TO BE CALCULATED ACROSS +0, +1, +2, NOT ACROSS TIME! TIME IS CONSIDERED VIA INTENSITY WEIGHTING
						float delta=auxScores[scanIndex][ppmIndex]-averageAuxScores[ppmIndex];
						precursorPPMVariance+=delta*delta;
					}
					if (auxScores.length>1) {
						precursorPPMVariance=precursorPPMVariance/(auxScores.length-1);
					}

					// rank in sortedBGScores
					int indexInSortedBGScores=Arrays.binarySearch(sortedBGScores, goodStripes.get(i).x);
					if (indexInSortedBGScores<0) {
						indexInSortedBGScores=-(indexInSortedBGScores+1);
					}
					int rank=sortedBGScores.length-indexInSortedBGScores;

					// averaging forward, so current scan is actually the median
					// for half a window back
					int medianIndex=index+scanAveragingHalfWindow;
					FragmentScan medianStripe=stripes.get(medianIndex);
					float[] completeAuxArray=General.concatenate(new float[] {numAboveThresholdMatches[index], numMatches[index], midTime[index]}, averageAuxScores,
							new float[] {fragmentDeltaMassAverage, fragmentDeltaMassVariance, duration, maxIDP, midIDP, precursorPPMVariance, bgsubScores[index],
									sumZScores[index]/scanAveragingWindow, rank, rawScores[index]});
					
					result.addStripe(goodStripes.get(i).x/scanAveragingWindow, completeAuxArray, 0.0f, 0.0f, medianStripe);

					if (identifiedPeaks>getPecanSearchParameters().getNumberOfReportedPeaks()) {
						// keep N+1 peaks
						break;
					}
					identifiedPeaks++;

					// add to takenScans to make sure we don't get a shoulder of
					// this peak
					int lowerWindow=medianIndex-2*scanAveragingWindow; // can't pick
																	// anything
																	// in twice
																	// the peak
																	// width
					int upperWindow=medianIndex+2*scanAveragingWindow; // +1 to
																	// account
																	// for the
																	// window
																	// boundary
					for (int j=lowerWindow; j<upperWindow; j++) {
						takenScans.add(j);
					}
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

			result.setTrace(new SortLaterXYTrace(scoreMap, GraphType.line, entry.getPeptideModSeq()));

			if (result.size()>0&&resultsQueue!=null) {
				try {
					resultsQueue.put(result);
				} catch (InterruptedException ie) {
					Logger.errorLine("Analysis interrupted!");
					Logger.errorException(ie);
					return Nothing.NOTHING;
				}
			}
		}
		return Nothing.NOTHING;
	}
}
