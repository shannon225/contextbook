package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;

import com.sun.xml.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.io.LineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.LineParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.math.FloatPair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RunningMedianWarper;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Sigmoid;
import edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest.RocPlot;
import gnu.trove.list.array.TFloatArrayList;

public class MProphet implements Runnable {
	private static final String DELIM = "\t";
	private final float peptideFDRThreshold;
	private final MProphetExecutionData settings;
	private final AminoAcidConstants aaConstants;
	
	private Throwable error;
	private Pair<ArrayList<PercolatorPeptide>, Float> result;

	public MProphet(MProphetExecutionData settings, float peptideFDRThreshold, AminoAcidConstants aaConstants) {
		this.settings = settings;
		this.peptideFDRThreshold=peptideFDRThreshold;
		this.aaConstants=aaConstants;
	}
	
	public static Pair<ArrayList<PercolatorPeptide>, Float> executeMProphetTSV(MProphetExecutionData commandData, float threshold, AminoAcidConstants aaConstants, int round) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		MProphet prophet=new MProphet(commandData, threshold, aaConstants);
		prophet.run();
		return prophet.getPeptides();
	}

	@Override
	public void run() {
		File file=settings.getInputTSV();

		try {
			MProphetDataset data = parseFeatureFile(file);
			result = calculateProbabilities(data);

		} catch (Throwable t) {
			Logger.errorLine("Error performing mProphet!");
			Logger.errorException(t);

			this.error = t;
		}
	}

	private Pair<ArrayList<PercolatorPeptide>, Float> calculateProbabilities(MProphetDataset dataset) throws EncyclopediaException {
//		System.out.println("Features: "+dataset.featureNames.size());
//		System.out.println("Targets: "+dataset.targetPeptideData.size());
//		System.out.println("Decoys: "+dataset.decoyPeptideData.size());
		LinearDiscriminantAnalysis lda=LinearDiscriminantAnalysis.buildModel(dataset.getTargetData(), dataset.getDecoyData());

//		double[] coefficients=lda.getCoefficients();
//		for (int i = 0; i < coefficients.length; i++) {
//			System.out.println(dataset.featureNames.get(i)+":\t"+coefficients[i]);
//		}
//		System.out.println("c:\t"+lda.getConstant());

		ArrayList<ScoredObject<MProphetData>> scoredDataset=new ArrayList<>();
		TFloatArrayList targetScores=new TFloatArrayList();
		TFloatArrayList decoyScores=new TFloatArrayList();
		for (MProphetData data : dataset.allData()) {
			float score=lda.getScore(data.data);
			scoredDataset.add(new ScoredObject<MProphet.MProphetData>(score, data));
			if (data.isDecoy) {
				decoyScores.add(score);
			} else {
				targetScores.add(score);
			}
		}

		// 1 stdev below the mean (still in the increasing section) or 10000 values, which ever is bigger 
		float b = Math.min(10000f/decoyScores.size(), 0.309f); // max at 0.5 stdevs below mean (if less than 30k peptides)
		b=Math.max(b, Math.min(1000f/decoyScores.size(), 0.5f)); // max at median (if less than 2000 peptides)
		float medianDecoy=QuickMedian.select(decoyScores.toArray(), Math.max(0.1586f, b));
		
		int bestBinCount=0;
		float pi0Estimate=0.0f;
		float bestSumSquaredErrors=Float.MAX_VALUE;
		for (int binCount = 10; binCount < 100; binCount++) {
			FloatPair weightedAverageData = getPi0Estimate(targetScores, decoyScores, medianDecoy, binCount);
			//System.out.println(binCount+"\t"+weightedAverageData.getOne()+"\t"+weightedAverageData.getTwo());
			
			if (weightedAverageData.getTwo()<bestSumSquaredErrors) {
				bestSumSquaredErrors=weightedAverageData.getTwo();
				pi0Estimate=weightedAverageData.getOne();
				bestBinCount=binCount;
			}
		}
		float pi0Prob = pi0Estimate/0.5f;

		ArrayList<XYPoint>[] points=PivotTableGenerator.createPivotTables(new float[][] {targetScores.toArray(), decoyScores.toArray()}, true);
		ArrayList<XYPoint> targets = points[0];
		ArrayList<XYPoint> decoys = points[1];

		ArrayList<XYPoint> scaledDecoys=new ArrayList<XYPoint>();
		ArrayList<XYPoint> delta=new ArrayList<XYPoint>();
		ArrayList<XYPoint> thresholdedDelta=new ArrayList<XYPoint>();
		float maxHistogram=0.0f;
		for (int i = 0; i < targets.size(); i++) {
			XYPoint target = targets.get(i);
			XYPoint decoy = decoys.get(i);
			maxHistogram=(float)Math.max(maxHistogram, target.y);
			double deltaRatio=target.y/(target.y+decoy.y);
			delta.add(new XYPoint(target.x, deltaRatio));
			thresholdedDelta.add(new XYPoint(target.x, Math.max(pi0Estimate, deltaRatio)));
			scaledDecoys.add(new XYPoint(decoy.x, decoy.y*pi0Prob));
		}
		int order=Math.max(3, Math.round(delta.size()/50f));
		RunningMedianWarper pepValueFunction=new RunningMedianWarper(thresholdedDelta, order, true);
		
		// calculate FDRs using T/D
		Collections.sort(scoredDataset);
		Collections.reverse(scoredDataset);
		float targetCount=0f;
		float decoyCount=0f;
		ArrayList<XYPoint> fdrCalc=new ArrayList<XYPoint>();
		for (ScoredObject<MProphetData> scoredData : scoredDataset) {
			float score=scoredData.getScore();
			boolean isDecoy=scoredData.y.isDecoy;
			
			if (isDecoy) {
				decoyCount+=pi0Prob;
			} else {
				targetCount++;
			}
			fdrCalc.add(new XYPoint(score, decoyCount/targetCount));
		}

		// calculate Q-values from FDRs
		ArrayList<XYPoint> qValueCalc=new ArrayList<XYPoint>();
		float movingQ=1.0f;
		for (int i = fdrCalc.size()-1; i >=0; i--) {
			XYPoint xyPoint = fdrCalc.get(i);
			movingQ=Math.min(movingQ, (float)xyPoint.y);
			qValueCalc.add(new XYPoint(xyPoint.x, movingQ));
		}
		LinearInterpolatedFunction qValueFunc=new LinearInterpolatedFunction(qValueCalc);

		// Find target peptides, estimate PEPs, and write files
		ArrayList<PercolatorPeptide> detectedPeptides=new ArrayList<>();
		float minScore=Float.MAX_VALUE;
		try {
			PrintWriter targetWriter=new PrintWriter(settings.getPeptideOutputFile(), "UTF-8");
			PrintWriter decoyWriter=new PrintWriter(settings.getPeptideDecoyFile(), "UTF-8");
			
			targetWriter.println("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds");
			decoyWriter.println("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds");
			
			for (ScoredObject<MProphetData> scoredData : scoredDataset) {
				float score=scoredData.getScore();
				
				float qValue=qValueFunc.getYValue(score);
				float posteriorErrorProb=pepValueFunction.getYValue(score);
				if (qValue<=peptideFDRThreshold&&!scoredData.y.isDecoy) {
					if (score<minScore) {
						minScore=score;
					}
					PercolatorPeptide pep=new PercolatorPeptide(scoredData.y.id, scoredData.y.protein, qValue, posteriorErrorProb, aaConstants);
					detectedPeptides.add(pep);
				}
				
				if (scoredData.y.isDecoy) {
					decoyWriter.println(scoredData.y.id+"\t"+score+"\t"+qValue+"\t"+posteriorErrorProb+"\t"+"-."+scoredData.y.sequence+".-"+"\t"+scoredData.y.protein);
				} else {
					targetWriter.println(scoredData.y.id+"\t"+score+"\t"+qValue+"\t"+posteriorErrorProb+"\t"+"-."+scoredData.y.sequence+".-"+"\t"+scoredData.y.protein);
				}
			}
			targetWriter.println("pi_0="+pi0Prob);
			decoyWriter.println("pi_0="+pi0Prob);

			targetWriter.flush();
			decoyWriter.flush();
			targetWriter.close();
			decoyWriter.close();
			
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: " + settings.getPeptideOutputFile().getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: " + settings.getPeptideOutputFile().getAbsolutePath(), e);
		}

//		XYTrace scoreMaxThresholdTrace=new XYTrace(new double[] {minScore, minScore}, new double[] {0, maxHistogram}, GraphType.dashedline, "Threshold");
//		XYTrace scoreThresholdTrace=new XYTrace(new double[] {minScore, minScore}, new double[] {0, 1}, GraphType.dashedline, "Threshold");
//		
//		XYTrace targetTrace=new XYTrace(targets, GraphType.line, "Target");
//		XYTrace decoyTrace=new XYTrace(scaledDecoys, GraphType.line, "Decoy*pi0");
//		Charter.launchChart("LDA Score", "Count", true, targetTrace, decoyTrace, scoreMaxThresholdTrace);
//
//		XYTrace fdrTrace = new XYTrace(fdrCalc, GraphType.line, "FDR");
//		XYTrace qvalueTrace=new XYTrace(qValueCalc, GraphType.dashedline, "Q-Value");
//		Charter.launchChart("LDA Score", "FDR/Q-Value", true, fdrTrace, qvalueTrace, scoreThresholdTrace);
//
//		XYTrace ratioTrace = new XYTrace(delta, GraphType.line, "Ratio");
//		XYTrace curveFit=new XYTrace(pepValueFunction.getKnots(), GraphType.dashedline, "Fit");
//		XYTrace piZeroTrace=new XYTrace(new double[] {targets.get(0).x, medianDecoy}, new double[] {pi0Estimate, pi0Estimate}, GraphType.dashedline, "Pi0");
//		Charter.launchChart("LDA Score", "Ratio", true, ratioTrace, curveFit, scoreThresholdTrace, piZeroTrace);
//		
//		try {Thread.sleep(1000000000);} catch (Exception e) {} // FIXME
		
		Pair<ArrayList<PercolatorPeptide>, Float> thisResult=new Pair<ArrayList<PercolatorPeptide>, Float>(detectedPeptides, pi0Prob);
		return thisResult;
	}

	private FloatPair getPi0Estimate(TFloatArrayList targetScores, TFloatArrayList decoyScores, float medianDecoy, int binCount) {
		ArrayList<XYPoint>[] localpoints=PivotTableGenerator.createPivotTables(new float[][] {targetScores.toArray(), decoyScores.toArray()}, true, binCount);
		ArrayList<XYPoint> localtargets = localpoints[0];
		ArrayList<XYPoint> localdecoys = localpoints[1];
		
		float sumWeights=0;
		float weightedAverage=0;
		ArrayList<XYPoint> localdelta=new ArrayList<XYPoint>();
		
		for (int i = 0; i < localtargets.size(); i++) {
			XYPoint target = localtargets.get(i);
			XYPoint decoy = localdecoys.get(i);
			double x=target.x;
			double deltaRatio=target.y/(target.y+decoy.y);
			localdelta.add(new XYPoint(x, deltaRatio));
			
			if (x<medianDecoy) {
				sumWeights+=decoy.y;
				weightedAverage+=decoy.y*deltaRatio;
			}
		}
		weightedAverage=weightedAverage/sumWeights;

		float sumSquaredErrors=0.0f;
		int n=0;
		for (int i = 0; i < localtargets.size(); i++) {
			XYPoint target = localtargets.get(i);
			double x=target.x;
			if (x<medianDecoy) {
				n++;
				XYPoint decoy = localdecoys.get(i);
				double deltaRatio=target.y/(target.y+decoy.y);
				double delta=weightedAverage-deltaRatio;
				sumSquaredErrors+=delta*delta*decoy.y;
			}
		}
		return new FloatPair(weightedAverage, sumSquaredErrors);
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
		return roc;
	}

	private MProphetDataset parseFeatureFile(File file) throws EncyclopediaException {
		String[] columnNames=null;
		int idIndex=0;
		int labelIndex=0;
		int scanNrIndex=0;
		int sequenceIndex=0;
		int proteinIndex=0;
		
		boolean[] isFeature=null;
		ArrayList<String> featureNames=new ArrayList<>();

		Logger.logLine("Parsing header for input file "+file.getName());
		try {
			BufferedReader in=new BufferedReader(new FileReader(file));
			String header=in.readLine();
			columnNames=header.split(DELIM, -1);
			isFeature=new boolean[columnNames.length];
			for (int i = 0; i < columnNames.length; i++) {
				if ("id".equals(columnNames[i])) {
					idIndex=i;
				} else if ("Label".equals(columnNames[i])) {
					labelIndex=i; // either TD or Label
				} else if ("TD".equals(columnNames[i])) {
					labelIndex=i; // either TD or Label
				} else if ("ScanNr".equals(columnNames[i])) {
					scanNrIndex=i;
				} else if ("sequence".equals(columnNames[i])) {
					sequenceIndex=i;
				} else if ("Proteins".equals(columnNames[i])) {
					proteinIndex=i;
				} else if ("pepLength".equals(columnNames[i])) {
					// skip
				} else if ("charge1".equals(columnNames[i])) {
					// skip
				} else if ("charge2".equals(columnNames[i])) {
					// skip
				} else if ("charge3".equals(columnNames[i])) {
					// skip
				} else if ("charge4".equals(columnNames[i])) {
					// skip
				} else if ("precursorMass".equals(columnNames[i])) {
					// skip
				} else if ("RTinMin".equals(columnNames[i])) {
					// skip
				} else if ("midTime".equals(columnNames[i])) {
					// skip
				} else {
					featureNames.add(columnNames[i]);
					isFeature[i]=true;
				}
			}
			in.close();
			Logger.logLine("Found indicies for "+featureNames.size()+" features: ["+General.toString(featureNames)+"]");

			final int idIndexFinal=idIndex;
			final int labelIndexFinal=labelIndex;
			final int sequenceIndexFinal=sequenceIndex;
			final int proteinIndexFinal=proteinIndex;
			final boolean[] isFeatureFinal=isFeature;
			final String[] columnNamesFinal=columnNames;
			
			ArrayList<MProphetData> peptideData=new ArrayList<>();
			LineParserMuscle muscle = new LineParserMuscle() {
				boolean isFirst=true;
				
				@Override
				public void processRow(String row) {
					if (isFirst) {
						// skip header
						isFirst=false;
						return;
					}
					String[] values=row.split(DELIM, -1);
					boolean isDecoy;
					
					try {
						isDecoy=Integer.parseInt(values[labelIndexFinal])<0;
					} catch (Exception e) {
						Logger.errorLine("Error parsing ["+values[labelIndexFinal]+"] as "+columnNamesFinal[labelIndexFinal]+" (index "+labelIndexFinal+")!");
						Logger.errorException(e);

						MProphet.this.error = e;
						throw e;
					}
					
					TFloatArrayList features=new TFloatArrayList();
					for (int j = 0; j < isFeatureFinal.length; j++) {
						if (isFeatureFinal[j]) {
							try {
								features.add(Float.parseFloat(values[j]));
							} catch (Exception e) {
								Logger.errorLine("Error parsing ["+values[j]+"] as "+columnNamesFinal[j]+" (index "+j+")!");
								Logger.errorException(e);

								MProphet.this.error = e;
								throw e;
							}
						}
					}
					peptideData.add(new MProphetData(values[idIndexFinal], values[sequenceIndexFinal], values[proteinIndexFinal], features.toArray(), isDecoy));
					
				}
				
				@Override
				public void cleanup() {
				}
			};
			LineParser.parseFile(settings.getInputTSV(), muscle);

			return new MProphetDataset(featureNames, peptideData);
			
		} catch (Throwable t) {
			Logger.errorLine("Error performing mProphet!");
			Logger.errorException(t);

			this.error = t;
			throw new EncyclopediaException(t);
		}
	}

	public boolean hadError() {
		return null != error;
	}

	public Throwable getError() {
		return error;
	}
	
	public Pair<ArrayList<PercolatorPeptide>, Float> getPeptides() {
		return result;
	}
	
	protected class MProphetDataset {
		private final ArrayList<String> featureNames;
		private final ArrayList<MProphetData> targetPeptideData;
		private final ArrayList<MProphetData> decoyPeptideData;
		public MProphetDataset(ArrayList<String> featureNames, ArrayList<MProphetData> peptideData) {
			this.featureNames = featureNames;
			targetPeptideData=new ArrayList<MProphet.MProphetData>();
			decoyPeptideData=new ArrayList<MProphet.MProphetData>();
			for (MProphetData mProphetData : peptideData) {
				if (mProphetData.isDecoy) {
					decoyPeptideData.add(mProphetData);
				} else {
					targetPeptideData.add(mProphetData);
				}
			}
		}
		
		public ArrayList<MProphetData> allData() {
			ArrayList<MProphetData> dataset=new ArrayList<>();
			dataset.addAll(targetPeptideData);
			dataset.addAll(decoyPeptideData);
			return dataset;
		}
		
		public ArrayList<float[]> getTargetData() {
			return getDataset(targetPeptideData);
		}
		
		public ArrayList<float[]> getDecoyData() {
			return getDataset(decoyPeptideData);
		}
		
		private ArrayList<float[]> getDataset(ArrayList<MProphetData> dataset) {
			ArrayList<float[]> data=new ArrayList<float[]>();
			for (MProphetData mProphetData : dataset) {
				data.add(mProphetData.data);
			}
			return data;
		}
		
		public ArrayList<MProphetData> getTargetPeptides() {
			return targetPeptideData;
		}
	}
	
	protected class MProphetData implements Comparable<MProphetData> {
		private final String id;
		private final String sequence;
		private final String protein;
		private final float[] data;
		private final boolean isDecoy;
		
		public MProphetData(String id, String sequence, String protein, float[] data, boolean isDecoy) {
			this.id = id;
			this.sequence = sequence;
			this.protein = protein;
			this.data = data;
			this.isDecoy=isDecoy;
		}
		
		@Override
		public int compareTo(MProphetData o) {
			return id.compareTo(o.id);
		}
		@Override
		public int hashCode() {
			return id.hashCode();
		}
		@Override
		public boolean equals(Object obj) {
			if (obj instanceof MProphetData) return compareTo((MProphetData)obj)==0;
			return false;
		}
	}
}
