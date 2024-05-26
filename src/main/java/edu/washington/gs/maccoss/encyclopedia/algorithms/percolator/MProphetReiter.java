package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

public class MProphetReiter implements Runnable {
	private final float peptideFDRThreshold;
	private final MProphetExecutionData settings;
	private final AminoAcidConstants aaConstants;
	
	private Throwable error;
	private MProphetResult result;
	
	public MProphetReiter(MProphetExecutionData settings, float peptideFDRThreshold, AminoAcidConstants aaConstants) {
		this.settings = settings;
		this.peptideFDRThreshold=peptideFDRThreshold;
		this.aaConstants=aaConstants;
	}
	
	public static MProphetResult executeMProphetTSV(MProphetExecutionData commandData, float threshold, AminoAcidConstants aaConstants, int round) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		MProphetReiter prophet=new MProphetReiter(commandData, threshold, aaConstants);
		prophet.run();
		return prophet.getResult();
	}

	@Override
	public void run() {
		File file=settings.getInputTSV();

		try {
			MProphetDataset data = LimitedMProphet.parseFeatureFile(file, settings);
			result = calculateProbabilities(data);

		} catch (Throwable t) {
			Logger.errorLine("Error performing mProphet!");
			Logger.errorException(t);

			this.error = t;
		}
	}

	private MProphetResult calculateProbabilities(MProphetDataset dataset) throws EncyclopediaException {
		int randomSeed=RandomGenerator.randomInt(1);
		int iterationCount = 50;
		int numIterationsPerCalculation = 10;
		
		ArrayList<LinearDiscriminantAnalysis> models=new ArrayList<LinearDiscriminantAnalysis>();
		for (int n = 0; n < iterationCount; n++) {
			randomSeed=RandomGenerator.randomInt(randomSeed);
			MProphetDataset[] folds=MProphetDataset.splitKFold(dataset, 2, randomSeed);
			
			MProphetDataset trainingDataset=folds[0];
			MProphetDataset testingDataset=folds[1];
			
			LinearDiscriminantAnalysis lda=null;
			int best=0;
			for (int i = 0; i < numIterationsPerCalculation; i++) {
				float targetFDR=0.01f;
				if (i==0) {
					targetFDR=0.15f;
				}
				ArrayList<ScoredMProphetData> data=trainingDataset.getPassingTargets(Optional.ofNullable(lda), targetFDR).x;
				
				if (data.size()<best) {
					break;
				}
				best=data.size();
				lda=LinearDiscriminantAnalysis.buildModel(MProphetDataset.getScoredData(data), trainingDataset.getDecoyData());
			}
			
			models.add(lda);
			Pair<ArrayList<ScoredMProphetData>, Float> data=testingDataset.getPassingTargets(Optional.ofNullable(lda), 0.01f);
			System.out.println("Iteration "+(n+1)+": "+data.x.size()+"/"+testingDataset.getTargetData().size()+" passing, pi0:"+data.y);
		}
		
		LinearDiscriminantAnalysis averageModel=LinearDiscriminantAnalysis.average(models);

		Pair<ArrayList<ScoredMProphetData>, Float> finalData=dataset.getPassingTargets(Optional.ofNullable(averageModel), Float.MAX_VALUE);
		int passingCount=0;
		for (ScoredMProphetData data : finalData.x) {
			if (data.fdr<0.01) passingCount++;
		}
		
		System.out.println("Final model: "+passingCount+"/"+dataset.getTargetData().size()+" passing, pi0:"+finalData.y);
		Pair<ArrayList<ScoredMProphetData>, Float> finalDecoyData=dataset.getPassingTargets(Optional.ofNullable(averageModel), Float.MAX_VALUE, true);
		
		ArrayList<ScoredMProphetData> allData=new ArrayList<ScoredMProphetData>();
		allData.addAll(finalData.x);
		allData.addAll(finalDecoyData.x);
			
		ArrayList<PercolatorPeptide> detectedPeptides=new ArrayList<>();
		float minScore=Float.MAX_VALUE;
		try {
			PrintWriter targetWriter=new PrintWriter(settings.getPeptideOutputFile(), "UTF-8");
			PrintWriter decoyWriter=new PrintWriter(settings.getPeptideDecoyFile(), "UTF-8");
			
			targetWriter.println("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds");
			decoyWriter.println("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds");
			
			for (ScoredMProphetData scoredData : allData) {
				float score=scoredData.getScore();
				
				float qValue=(float)scoredData.getFDR();
				float posteriorErrorProb=(float)scoredData.getLocalFDR();
				if (qValue<=peptideFDRThreshold&&!scoredData.getData().isDecoy()) {
					if (score<minScore) {
						minScore=score;
					}
					PercolatorPeptide pep=new PercolatorPeptide(scoredData.getData().getId(), scoredData.getData().getProtein(), qValue, posteriorErrorProb, aaConstants);
					detectedPeptides.add(pep);
				}
				
				if (scoredData.getData().isDecoy()) {
					decoyWriter.println(scoredData.getData().getId()+"\t"+score+"\t"+qValue+"\t"+posteriorErrorProb+"\t"+"-."+scoredData.getData().getSequence()+".-"+"\t"+scoredData.getData().getProtein());
				} else {
					targetWriter.println(scoredData.getData().getId()+"\t"+score+"\t"+qValue+"\t"+posteriorErrorProb+"\t"+"-."+scoredData.getData().getSequence()+".-"+"\t"+scoredData.getData().getProtein());
				}
			}
			targetWriter.println("pi_0="+finalData.y);
			decoyWriter.println("pi_0="+finalData.y);

			targetWriter.flush();
			decoyWriter.flush();
			targetWriter.close();
			decoyWriter.close();
			
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: " + settings.getPeptideOutputFile().getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: " + settings.getPeptideOutputFile().getAbsolutePath(), e);
		}
		
		return new MProphetResult(detectedPeptides, averageModel, dataset.getFeatureNames(), finalData.y);
	}

	public boolean hadError() {
		return null != error;
	}

	public Throwable getError() {
		return error;
	}
	
	public MProphetResult getResult() {
		return result;
	}
}
