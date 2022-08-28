package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
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
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest.RocPlot;
import gnu.trove.list.array.TFloatArrayList;

public class MProphet implements Runnable {
	private static final String DELIM = "\t";
	private final float peptideFDRThreshold;
	private final MProphetExecutionData settings;
	private Throwable error;
	private Pair<ArrayList<PercolatorPeptide>, Float> result;

	public MProphet(MProphetExecutionData settings, float peptideFDRThreshold) {
		this.settings = settings;
		this.peptideFDRThreshold=peptideFDRThreshold;
	}
	
	public static Pair<ArrayList<PercolatorPeptide>, Float> executeMProphetTSV(MProphetExecutionData commandData, float threshold, AminoAcidConstants aaConstants, int round) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		MProphet prophet=new MProphet(commandData, threshold);
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
		System.out.println("Features: "+dataset.featureNames.size());
		System.out.println("Targets: "+dataset.targetPeptideData.size());
		System.out.println("Decoys: "+dataset.decoyPeptideData.size());
		LinearDiscriminantAnalysis lda=LinearDiscriminantAnalysis.buildModel(dataset.getTargetData(), dataset.getDecoyData());

		double[] coefficients=lda.getCoefficients();
		for (int i = 0; i < coefficients.length; i++) {
			System.out.println(dataset.featureNames.get(i)+":\t"+coefficients[i]);
		}
		System.out.println("c:\t"+lda.getConstant());

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

		ArrayList<XYPoint>[] points=PivotTableGenerator.createPivotTables(new float[][] {targetScores.toArray(), decoyScores.toArray()}, true);
		
		float medianDecoy=QuickMedian.median(decoyScores.toArray());
		
		float sumWeights=0;
		float weightedAverage=0;
		ArrayList<XYPoint> delta=new ArrayList<XYPoint>();
		for (int i = 0; i < points[0].size(); i++) {
			XYPoint target = points[0].get(i);
			XYPoint decoy = points[1].get(i);
			double x=target.x;
			double deltaRatio=target.y/(target.y+decoy.y);
			delta.add(new XYPoint(x, deltaRatio));
			
			if (x<medianDecoy) {
				sumWeights+=decoy.y;
				weightedAverage+=decoy.y*deltaRatio;
			}
		}
		weightedAverage=weightedAverage/sumWeights;
		
		XYTraceInterface[] traces=new XYTraceInterface[2];
		traces[0]=new XYTrace(points[0], GraphType.line, "Target");
		traces[1]=new XYTrace(points[1], GraphType.line, "Decoy");
		Charter.launchChart("LDA Score", "Count", true, traces);

		XYTrace ratioTrace = new XYTrace(delta, GraphType.line, "Ratio");
		XYTrace piZeroTrace=new XYTrace(new double[] {points[0].get(0).x, medianDecoy}, new double[] {weightedAverage, weightedAverage}, GraphType.dashedline, "Pi0");
		Charter.launchChart("LDA Score", "Ratio", true, ratioTrace, piZeroTrace);
		
		try {
		Thread.sleep(1000000000);
		} catch (Exception e) {}
		
		Pair<ArrayList<PercolatorPeptide>, Float> thisResult=new Pair<ArrayList<PercolatorPeptide>, Float>(null, null);
		return thisResult;
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
		throw new RuntimeException("NOT IMPLEMENTED YET");
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
