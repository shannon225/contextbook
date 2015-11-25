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
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.TFloatFloatHashMap;
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
			int requiredNumAboveThreshold=(int)(0.5f*entry.getPeptideSeq().length());
			
			
			float[] rawScores=new float[super.stripes.size()];
			float[] bgsubScores=new float[super.stripes.size()];
			float[][] fragmentTraces=new float[entry.getMassArray().length][];
			for (int i=0; i<fragmentTraces.length; i++) {
				fragmentTraces[i]=new float[super.stripes.size()];
			}
			
			for (int i=0; i<super.stripes.size(); i++) {
				Stripe stripe=super.stripes.get(i);
				float[] scores=getScorer().getIndividualPeakScores(entry, stripe, false);
				for (int j=0; j<scores.length; j++) {
					fragmentTraces[j][i]=scores[j];
				}
				rawScores[i]+=getScorer().score(entry, stripe);
				
				float rt=stripe.getScanStartTime();
				XYPoint meanStdev=background.get((double)rt);
				if (meanStdev!=null) {
					bgsubScores[i]=(float)(rawScores[i]-meanStdev.x);
				} else {
					bgsubScores[i]=rawScores[i];
				}
			}

			int scanAveragingWindow=2*scanAveragingMargin+1;
			int scanExcludingWindow=2*scanAveragingWindow+1;

			float[] sumRawScores=movingSum(rawScores, scanAveragingWindow);
			float[] sumBgsubScores=movingSum(bgsubScores, scanAveragingWindow);
			float[][] sumFragmentTraces=new float[entry.getIntensityArray().length][];
			for (int i=0; i<sumFragmentTraces.length; i++) {
				sumFragmentTraces[i]=movingSum(fragmentTraces[i], scanAveragingWindow);
			}

			ArrayList<ScoredObject<IndexedObject<Stripe>>> goodStripes=new ArrayList<ScoredObject<IndexedObject<Stripe>>>();
			int[] numAboveThresholdMatches=new int[sumRawScores.length];
			int[] numMatches=new int[sumRawScores.length];
			for (int i=0; i<numAboveThresholdMatches.length; i++) {
				// NOTE: this seems questionable that unnormalized intensities are used for individual scores while normalized intensities are used for total scores. -BCS
				float threshold=sumRawScores[i]/(entry.getPeptideSeq().length()+1);
				for (int j=0; j<sumFragmentTraces.length; j++) {
					if (sumFragmentTraces[j][i]>=threshold) {
						numAboveThresholdMatches[i]++;
					}
					if (sumFragmentTraces[j][i]>=threshold) {
						numMatches[i]++;
					}
				}
				
				/*if (super.stripes.get(i).getScanStartTime()>60f*48.85&&super.stripes.get(i).getScanStartTime()<60f*49.20f) {
					System.out.println(entry.getPeptideSeq()+"\t"+(super.stripes.get(i).getScanStartTime()/60f)+"\t"+numAboveThresholdMatches[i]+"/"+requiredNumAboveThreshold+"\t"+sumBgsubScores[i]);
					for (int j=0; j<sumFragmentTraces.length; j++) {
						System.out.println("\t"+sumFragmentTraces[j][i]+"\t>= "+threshold);
					}
				}*/
				/*
				if (super.stripes.get(i).getScanStartTime()>60f*42.32f&&super.stripes.get(i).getScanStartTime()<60f*42.33f) {
					System.out.println(entry.getPeptideSeq()+"\t"+(super.stripes.get(i).getScanStartTime()/60f)+"\t"+numAboveThresholdMatches[i]+"/"+requiredNumAboveThreshold+"\t"+sumBgsubScores[i]);
					for (int j=0; j<sumFragmentTraces.length; j++) {
						System.out.println("\t"+sumFragmentTraces[j][i]+"\t>= "+threshold);
					}
				}*/
				
				if (numAboveThresholdMatches[i]>requiredNumAboveThreshold) {
					goodStripes.add(new ScoredObject<IndexedObject<Stripe>>(sumBgsubScores[i], new IndexedObject<Stripe>(i, stripes.get(i))));
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
					for (int j=0; j<scanExcludingWindow; j++) {
						int index=stripe.x-scanAveragingWindow+j;
						if (index>=0&&index<sumBgsubScores.length) {
							takenScans.add(index);
						}
					}
					for (int j=0; j<scanAveragingWindow; j++) {
						int index=stripe.x-scanAveragingMargin+j;
						if (index>=0&&index<stripes.size()) {
							float[] auxScores=getScorer().auxScore(entry, stripes.get(index));
							if (averageAuxScores==null) {
								averageAuxScores=auxScores;
							} else {
								averageAuxScores=General.add(averageAuxScores, auxScores);
							}
						}
					}
					averageAuxScores=General.multiply(averageAuxScores, 1.0f/scanAveragingWindow);
					
					result.addStripe(goodStripes.get(i).x, General.concatenate(averageAuxScores, numAboveThresholdMatches[stripe.x], numMatches[stripe.x], sumRawScores[stripe.x]), stripe.y);
				}
			}
			
			TFloatFloatHashMap scoreMap=new TFloatFloatHashMap();
			for (int i=0; i<super.stripes.size(); i++) {
				if (numAboveThresholdMatches[i]>=requiredNumAboveThreshold) {
					scoreMap.put(super.stripes.get(i).getScanStartTime(), sumBgsubScores[i]);
				} else {
					scoreMap.put(super.stripes.get(i).getScanStartTime(), 0.0f);
				}
			}
			
			result.setTrace(new XYTrace(scoreMap, GraphType.line, entry.getPeptideModSeq()));
			map.put(entry, result);
		}
		return map;
	}
	
	public static float[] movingSum(float[] scores, int scanAveragingWindow) {
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
