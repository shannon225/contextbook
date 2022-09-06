package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;

import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest.RocPlot;

public class LinearDiscriminantAnalysisTest {
	
	
	public static void main(String[] args) {
		File f=new File("/Users/searleb/Documents/damien/dda_library_search/scoring_tests/multiple_evalues/23aug2017_hela_serum_timecourse_pool_dda_001.dia.scribe.features.txt");
//		final String[] headers = new String[] { "e2Score", "evalue", "correlationToGaussian", "correlationToPrecursor",
//				"correlationToPlusOne", "isIntegratedSignal", "isIntegratedPrecursor", "numPeaksWithGoodCorrelation",
//				"numPeaksWithGreatCorrelation", "primary", "xCorrLib", "xCorrModel", "LogDotProduct",
//				"logWeightedDotProduct", "sumOfSquaredErrors", "weightedSumOfSquaredErrors", "numberOfMatchingPeaks",
//				"numberOfMatchingPeaksAboveThreshold", "averageAbsFragmentDeltaMass", "averageFragmentDeltaMasses",
//				"isotopeDotProduct", "averageAbsParentDeltaMass", "averageParentDeltaMass", "deltaRT",
//				"numMissedCleavage" };
		
		//id	Label	ScanNr	DotProduct	contrastAngle	logit	primary	xCorrLib	xCorrModel	lnInvSumOfSquaredErrors	numberOfMatchingPeaks	averageAbsFragmentDeltaMass	averageFragmentDeltaMasses	isotopeDotProduct	averageAbsParentDeltaMass	averageParentDeltaMass	percentBlankOverMono	numberPrecursorMatch	Sp	maxLadderLength	eValueXCorr	eValueTandem	eValueSumSquares	oneHotXCorr	oneHotTandem	oneHotSumSquares	numConsidered	deltaRT	numMissedCleavage	pepLength	charge1	charge2	charge3	charge4	precursorMz	precursorMass	RTinMin	sequence	Proteins

		
		final String[] headers = new String[] { "xCorrLib" , "primary", "lnInvSumOfSquaredErrors"};
		
//		File f=new File("/Users/searleb/Documents/teaching/encyclopedia/test/old_encyclopedia/23aug2017_hela_serum_timecourse_wide_1a.dia.features.txt");		
//		final String[] headers = new String[] { "primary", "xCorrLib", "xCorrModel", "LogDotProduct",
//				"logWeightedDotProduct", "sumOfSquaredErrors", "weightedSumOfSquaredErrors", "numberOfMatchingPeaks",
//				"numberOfMatchingPeaksAboveThreshold", "averageAbsFragmentDeltaMass", "averageFragmentDeltaMasses",
//				"isotopeDotProduct", "averageAbsParentDeltaMass", "averageParentDeltaMass", "eValue", "deltaRT",
//				"numMissedCleavage" };
		
		ArrayList<float[]> target=new ArrayList<>();
		ArrayList<float[]> decoy=new ArrayList<>();
		
		//final float[] sums=new float[headers.length];
		
		TableParser.parseTSV(f, new TableParserMuscle() {
			
			@Override
			public void processRow(Map<String, String> row) {
				int td=Integer.parseInt(row.get("Label"));
				
				float[] data=new float[headers.length];
				for (int i = 0; i < headers.length; i++) {
					data[i]=Float.parseFloat(row.get(headers[i]));
					//sums[i]+=Math.abs(data[i]);
				}
				
				if (td>0) {
					target.add(data);
				} else {
					decoy.add(data);
				}
			}
			
			@Override
			public void cleanup() {
				
			}
		});
		//int count=target.size()+decoy.size();
		
		if (true) {
			// calculate coefficients for current scores
			
			LinearDiscriminantAnalysis lda=LinearDiscriminantAnalysis.buildModel(target, decoy);
			double[] coefficients=lda.getCoefficients();
			for (int i = 0; i < coefficients.length; i++) {
				System.out.println(headers[i]+"\t"+coefficients[i]);
			}
			System.out.println("c: "+lda.getConstant());

//			System.out.println();
//			for (float[] data : target) {
//				System.out.println("1\t"+lda.getScore(data));
//			}
//			for (float[] data : decoy) {
//				System.out.println("-1\t"+lda.getScore(data));
//			}
			
//			0: -0.013720960802214612
//			1: -0.1075904065402189
//			2: 0.9523161866448213
//			c: -0.2812703285341387
			
		} else {
			// calculate masking to identify best features
		
			int[] lastMask=new int[headers.length];
			
			while (true) {
				float bestTPR=0;
				int bestCurrentMask=0;
				
				System.out.println("BEST CURRENT LIST:");
				for (int i = 0; i < lastMask.length; i++) {
					if (lastMask[i]>0) {
						System.out.print(" "+headers[i]);
					}
				}
				System.out.println();
				System.out.println();
				
				for (int m = 0; m < headers.length; m++) {
					if (lastMask[m]>0) continue;
					
					int[] currentMask=lastMask.clone();
					currentMask[m]=1;
					
					String mask=headers[m];	
					ArrayList<float[]> positive=mask(target, currentMask);
					ArrayList<float[]> negative=mask(decoy, currentMask);
					
					LinearDiscriminantAnalysis lda=LinearDiscriminantAnalysis.buildModel(positive, negative);
		//			double[] coefficients=lda.getCoefficients();
		//			for (int i = 0; i < coefficients.length; i++) {
		//				System.out.println(headers[i]+"\t"+coefficients[i]);
		//			}
		//			System.out.println("c: "+lda.getConstant());
					
					ArrayList<ScoredIndex> scores=new ArrayList<>();
					for (int j = 0; j < positive.size(); j++) {
						scores.add(new ScoredIndex(lda.getScore(positive.get(j)), 1));
					}
					for (int j = 0; j < negative.size(); j++) {
						scores.add(new ScoredIndex(lda.getScore(negative.get(j)), 0));
					}
			
					Collections.sort(scores);
					Collections.reverse(scores);
					
					RocPlot plot=getRocPlot(scores, positive.size(), negative.size());
					float tpr = plot.getTPR(0.01f)*positive.size();
					System.out.println("AUC: "+plot.getAUC(0.01f)+", numTargets: "+tpr+", eval cases: "+scores.size()+" --> "+mask);
					if (tpr>bestTPR) {
						bestTPR=tpr;
						bestCurrentMask=m;
					}
				}
				
				lastMask[bestCurrentMask]=1;
				
				if (General.sum(lastMask)==headers.length) break;
			}

		}
	}
	
	private static ArrayList<float[]> mask(ArrayList<float[]> data, int[] mask) {
		int maskLength=General.sum(mask);
		ArrayList<float[]> r=new ArrayList<float[]>();
		for (float[] fs : data) {
			float[] a=new float[maskLength];
			int index=0;
			for (int i = 0; i < fs.length; i++) {
				if (mask[i]>0) {
					a[index]=fs[i];
					index++;
				}
			}
			r.add(a);
		}
		return r;
	}
	
	private static RocPlot getRocPlot(ArrayList<ScoredIndex> scores, int totalPositives, int totalNegatives) {
		RocPlot roc = new RocPlot();
		int falsePositives = 0;
		int truePositives = 0;
		for (ScoredIndex r : scores) {
			if (r.y==1) {
				truePositives++;
			} else {
				falsePositives++;
			}

			float falsePositiveRate = falsePositives / (float) totalNegatives;
			float truePositiveRate = truePositives / (float) totalPositives;
			roc.addData(falsePositiveRate, truePositiveRate);
		}
		//System.out.println("ROC:\n"+roc.toString());
		return roc;
	}
}
