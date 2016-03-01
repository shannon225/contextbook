package edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest;




import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedianDouble;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.TFloatHashSet;

public class RandomForest {
	public static final float log2=(float)Math.log(2.0f);
	public static final int hardMinInEachGroup=9;

	public static RocPlot classify(ArrayList<Record> records, TIntIntHashMap enrichedIndicies) {
		
		// need to maintain the distribution of samples
		ArrayList<Record> positive=new ArrayList<Record>();
		ArrayList<Record> negative=new ArrayList<Record>();
		for (Record record : records) {
			if (record.state()) {
				positive.add(record);
			} else {
				negative.add(record);
			}
		}
		ArrayList<Record> testing=new ArrayList<Record>();
		float testFraction=0.2f;
		int positiveTestCount = (int)Math.max(hardMinInEachGroup, positive.size()*testFraction);
		positiveTestCount+=(int)(positiveTestCount/5.0f)*(float)Math.random(); // ranges from 10+/-1 (or 20%)
		for (int i=0; i<positiveTestCount; i++) {
			Record removed = positive.remove((int)(Math.random()*positive.size()));
			records.remove(removed);
			testing.add(removed);
		}
		int negativeTestCount = Math.min(Math.round((20.0f*positiveTestCount)/3.0f), (int)Math.max(hardMinInEachGroup, negative.size()*testFraction));
		negativeTestCount+=(int)(negativeTestCount/5.0f)*(float)Math.random(); // ranges from 10+/-1 (or 20%)
		for (int i=0; i<negativeTestCount; i++) {
			Record removed = negative.remove((int)(Math.random()*negative.size()));
			records.remove(removed);
			testing.add(removed);
		}
		//System.out.println(testFraction+"\t"+positiveTestCount+"\t"+negativeTestCount+"\t"+records.size()+"\t"+testing.size());
		
		return classify(records, testing, enrichedIndicies);
	}

	public static RocPlot classify(ArrayList<Record> records,
			ArrayList<Record> testing, TIntIntHashMap enrichedIndicies) {
		boolean[] useData=trimDataSet(records);
		trimDataSet(testing, useData);
		int useDataCount=0;
		for (int i = 0; i < useData.length; i++) {
			if (useData[i]) useDataCount++;
		}

		ArrayList<ScoredObject<TreeNode>> scoredTrees=new ArrayList<ScoredObject<TreeNode>>();
		Forest forest=new Forest();
		for (int i=1; i<=100; i++) {
			ArrayList<Record> bootstrap=new ArrayList<Record>();
			for (int j=0; j<records.size(); j++) {
				bootstrap.add(records.get((int)(Math.random()*records.size())));
			}

			ArrayList<Record> reserved=new ArrayList<Record>();
			for (Record record : records) {
				boolean used=false;
				for (Record inBootstrap : bootstrap) {
					if (record==inBootstrap) {
						used=true;
						break;
					}
				}
				if (!used) {
					reserved.add(record);
				}
			}

			TreeNode root=createNode(null, bootstrap);
			float probTreeReserved = getTreeAccuracy(root, reserved);
			for (int j = 0; j < 1000; j++) {
				TreeNode leaf=root.getRandomLeaf();
				TreeNode parent=leaf.getParent();
				if (parent==null) continue;

				boolean wasLeft=true;
				if (parent.getLeft() == leaf) {
					parent.setLeft(null);
				} else if (parent.getRight() == leaf) {
					parent.setRight(null);
					wasLeft=false;
				}

				float trimmedProbTreeReserved = getTreeAccuracy(root, reserved);
				if (trimmedProbTreeReserved>probTreeReserved) {
					// trimming helped
					probTreeReserved=trimmedProbTreeReserved;
				} else {
					// trimming didn't help, so reinstate
					if (wasLeft) {
						parent.setLeft(leaf);
					} else {
						parent.setRight(leaf);
					}
				}
			}
			
			scoredTrees.add(new ScoredObject<TreeNode>(probTreeReserved, root));
		}
		float averageReserved=0.0f;
		Collections.sort(scoredTrees);
		for (int i = 0; i < scoredTrees.size(); i++) {
			averageReserved+=scoredTrees.get(i).x;
			TreeNode root=scoredTrees.get(i).y;
			forest.add(root);

			/*
			int correctForest = 0;
			int incorrectForest = 0;

			for (Record record : testing) {
				float probability = forest.predictProbability(record);
				if (record.state() == probability > 0.5f) {
					correctForest++;
				} else {
					incorrectForest++;
				}
			}
			float probForest = correctForest / (float) (correctForest + incorrectForest);

			float probTree = getTreeAccuracy(root, testing);
			float probTraining = getTreeAccuracy(root, records);
			//System.out.println(i + "\t" + probTree + "\t" + probForest + "\t" + scoredTrees.get(i).score);
			 */
		}
		
		ArrayList<ScoredObject<Record>> scores=new ArrayList<ScoredObject<Record>>();
		int totalPositives=0;
		int totalNegatives=0;
		for (Record record : testing) {
			float probability=forest.predictProbability(record);
			scores.add(new ScoredObject<Record>(probability, record));
			if (record.state()) {
				totalPositives++;
			} else {
				totalNegatives++;
			}
		}

		Collections.sort(scores);
		Collections.reverse(scores);
		
		RocPlot roc = getRocPlot(scores, totalPositives, totalNegatives);
		
		// shuffle state and regenerate the roc
		float randomAUC=shuffleAndGetRocPlot(scores, totalPositives, totalNegatives).getAUC();
		
		forest.addDecisionIndicies(enrichedIndicies);
		
		// data collection about test
		int positiveTest=0;
		int negativeTest=0;
		for (Record test : testing) {
			if (test.state()) {
				positiveTest++;
			} else {
				negativeTest++;
			}
		}

		// data collection about training
		int positiveTrain=0;
		int negativeTrain=0;
		for (Record test : records) {
			if (test.state()) {
				positiveTrain++;
			} else {
				negativeTrain++;
			}
		}
		float target=positiveTrain/(float)(positiveTrain+negativeTrain);
		
		System.out.println(General.max(enrichedIndicies.values())+"\t"+positiveTest+"\t"+negativeTest+"\t"+positiveTrain+"\t"+negativeTrain+"\t"+target+"\t"+useDataCount+"\t"+(averageReserved/scoredTrees.size())+"\t"+roc.getAUC()+"\t"+randomAUC);
		return roc;
	}

	private static RocPlot shuffleAndGetRocPlot(ArrayList<ScoredObject<Record>> scores,
			int totalPositives, int totalNegatives) {
		ArrayList<ScoredObject<Record>> permuted=new ArrayList<ScoredObject<Record>>();
		Collections.shuffle(scores);
		ScoredObject<Record> prev=scores.get(scores.size()-1);
		for (int i = 0; i < scores.size(); i++) {
			ScoredObject<Record> orig=scores.get(i);
			ScoredObject<Record> record=new ScoredObject<Record>(orig.x, new Record(orig.y.getSampleName(), prev.y.state(), orig.y.getData()));
			permuted.add(record);
			prev=orig;
		}

		Collections.sort(permuted);
		Collections.reverse(permuted);
		
		RocPlot randomRoc = getRocPlot(permuted, totalPositives, totalNegatives);
		return randomRoc;
	}

	private static RocPlot getRocPlot(ArrayList<ScoredObject<Record>> scores,
			int totalPositives, int totalNegatives) {
		RocPlot roc=new RocPlot();
		int falsePositives=0;
		int truePositives=0;
		for (ScoredObject<Record> r : scores) {
			if (r.y.state()) {
				truePositives++;
			} else {
				falsePositives++;
			}
			
			float falsePositiveRate = falsePositives/(float)totalNegatives;
			float truePositiveRate = truePositives/(float)totalPositives;
			roc.addData(falsePositiveRate, truePositiveRate);
			//System.out.println(r.getR().state()+"\t"+r.getScore()+"\t"+falsePositives+"\t"+truePositives+"\t"+falsePositiveRate+"\t"+truePositiveRate);
		}
		return roc;
	}
	
	private static ArrayList<Record> trimDataSet(ArrayList<Record> records, boolean[] useData) {
		ArrayList<Record> newRecords=new ArrayList<Record>();
		for (Record record : records) {
			float[] data=record.getTrimmedData(useData);
			newRecords.add(new Record(record.getSampleName(), record.state(), data));
		}
		return newRecords;
	}
	
	private static boolean[] trimDataSet(ArrayList<Record> records) {
		MannWhitneyUTest test=new MannWhitneyUTest();
		
		double[] pValues=new double[records.get(0).getDataLength()];
		for (int g = 0; g < pValues.length; g++) {
			TDoubleArrayList s1=new TDoubleArrayList();
			TDoubleArrayList s2=new TDoubleArrayList();
			for (Record record : records) {
				if (record.state()) {
					s1.add(record.getDataValue(g));
				} else {
					s2.add(record.getDataValue(g));
				}
			}
			double[] ref = s1.toArray();
			double[] nref = s2.toArray();
			double pvalue=test.mannWhitneyUTest(ref, nref);
			pValues[g]=pvalue;
		}
		double target=QuickMedianDouble.select(pValues.clone(), 0.01);
		//System.out.println("TARGET: "+target+"\t"+sort[0]+"\t"+sort[sort.length-1]);
		boolean[] useData=new boolean[pValues.length];
		for (int g = 0; g < pValues.length; g++) {
			useData[g]=pValues[g]<=target;
		}
		
		trimDataSet(records, useData);
		return useData;
	}
	
	public static void testDataSetWithMannWhitney(ArrayList<Record> records) {
		MannWhitneyUTest test=new MannWhitneyUTest();
		
		double[] pValues=new double[records.get(0).getDataLength()];
		double[] ratios=new double[pValues.length];
		for (int g = 0; g < pValues.length; g++) {
			TDoubleArrayList s1=new TDoubleArrayList();
			TDoubleArrayList s2=new TDoubleArrayList();
			for (Record record : records) {
				if (record.state()) {
					s1.add(record.getDataValue(g));
				} else {
					s2.add(record.getDataValue(g));
				}
			}
			double[] ref = s1.toArray();
			double[] nref = s2.toArray();

			ratios[g]=Math.log10(General.mean(ref))-Math.log10(General.mean(nref));
			double pvalue;
			if (nref.length>ref.length*5) {
				TDoubleArrayList pvalues=new TDoubleArrayList();
				for (int i = 0; i < 100; i++) {
					TDoubleArrayList random=new TDoubleArrayList();
					for (int j = 0; j < ref.length; j++) {
						random.add(s2.get((int)Math.floor(Math.random()*s2.size())));
					}
					pvalues.add(test.mannWhitneyUTest(ref, random.toArray()));
				}
				pvalue=General.mean(pvalues.toArray());
			} else {
				pvalue=test.mannWhitneyUTest(ref, nref);
			}
			pValues[g]=pvalue;
			if (Double.isNaN(ratios[g])) ratios[g]=0.0;
			if (Double.isNaN(pValues[g])) pValues[g]=1.0;
		}
		double target=QuickMedianDouble.select(pValues.clone(), 0.01);
		//System.out.println("TARGET: "+target+"\t"+sort[0]+"\t"+sort[sort.length-1]);
		for (int g = 0; g < pValues.length; g++) {
			if (pValues[g]<=target) {
				System.out.println((g+1)+"\t"+(float)pValues[g]+"\t"+(float)ratios[g]);
			}
		}
		
	}

	private static float getTreeAccuracy(TreeNode root, ArrayList<Record> reserved) {
		int correctTree=0;
		int incorrectTree=0;
		for (Record record : reserved) {
			if (record.state()==root.predictValue(record)) {
				correctTree++;
			} else {
				incorrectTree++;
			}
		}
		float probTree=correctTree/(float)(correctTree+incorrectTree);
		return probTree;
	}
	
	public static TreeNode createNode(TreeNode parent, ArrayList<Record> records) {
		if (records.size()==0) return null;
		
		float prob=calculateProbability(records);
		if (Math.min(prob, 1.0f-prob)<0.1f) {
			// stop calculating when you're already biased 90% to one side
			return new TreeNode(parent, null, prob);
		} else if (records.size()>2&&Math.min(prob, 1.0f-prob)*records.size()<1f) {
			// also stop if you have between 3 and 10 records and your split is 1 to N-1
			return new TreeNode(parent, null, prob);
		}
		
		DecisionTreeSplitPoint point=chooseCutoff(records);
		if (point==null) {
			return new TreeNode(parent, null, prob);
		}
		ArrayList<Record>[] split=point.split(records);
		TreeNode node=new TreeNode(parent, point, prob);
		node.setLeft(createNode(node, split[0]));
		node.setRight(createNode(node, split[1]));
		return node;
	}

	public static DecisionTreeSplitPoint chooseCutoff(ArrayList<Record> records) {
		DecisionTreeSplitPoint bestDecision=null;
		
		int numberOfRandomFeatures=(int)Math.sqrt(records.get(0).getData().length);
		for (int f=0; f<numberOfRandomFeatures; f++) {
			int index=(int)(Math.random()*records.get(0).getData().length);
			
			float totalEntropy=calculateEntropy(records);
			TFloatHashSet scoreSet=new TFloatHashSet();
			for (Record record : records) {
				scoreSet.add(record.getData()[index]);
			}
			float[] scores=scoreSet.toArray();
			Arrays.sort(scores);
			// until scores.length-1 because we're using a gt threshold
			float cutoff=0.0f;
			float maxGain=-Float.MAX_VALUE;
			
			// try 100 random cutoffs, choose the best
			for (int i=0; i<100; i++) {
				float threshold=scores[(int)(scores.length*Math.random())];
				ArrayList<Record> left=new ArrayList<Record>();
				ArrayList<Record> right=new ArrayList<Record>();
				for (Record record : records) {
					if (record.getData()[index]>threshold) {
						right.add(record);
					} else {
						left.add(record);
					}
				}
				float leftEntropy=calculateEntropy(left);
				float rightEntropy=calculateEntropy(right);
				float gain=totalEntropy-(left.size()*leftEntropy/records.size()+right.size()*rightEntropy/records.size());
				if (gain>maxGain) {
					maxGain=gain;
					cutoff=threshold;
				}
			}
			DecisionTreeSplitPoint decision=new DecisionTreeSplitPoint(index, cutoff, maxGain);
			if (bestDecision==null) {
				bestDecision=decision;
			} else if (decision.gain>bestDecision.gain) {
				bestDecision=decision;
			}
		}
		return bestDecision;
	}

	public static float calculateEntropy(ArrayList<Record> records) {
		float pPos=calculateProbability(records);
		float pNeg=1.0f-pPos;
		return (-pPos*log2(pPos))-(pNeg*log2(pNeg));
	}

	private static float calculateProbability(ArrayList<Record> records) {
		int trues=0;
		int falses=0;
		for (Record record : records) {
			if (record.state()) {
				trues++;
			} else {
				falses++;
			}
		}
		float pPos=trues/(float)(trues+falses);
		return pPos;
	}

	public static float log2(float v) {
		if (v<=0.0f) return 0.0f;
		return (float)Math.log(v)/log2;
	}
}
