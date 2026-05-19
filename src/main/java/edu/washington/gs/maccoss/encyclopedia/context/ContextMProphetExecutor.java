package edu.washington.gs.maccoss.encyclopedia.context;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetDataset;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetFeatureReader;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetReiter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;

public class ContextMProphetExecutor {

	public static void main(String[] args) {
		
		// Map files 
		String outputPath = "C:/Users/m334793/Documents/Library/for_context_50perCycle/";
		String libraryPath = "C:/Users/m334793/Documents/Library/for_context_50perCycle/IL2_and_IL15_Combo.elib";
		String fastaPath = "C:/Users/m334793/Documents/Library/for_context_50perCycle/mus_musculus_reviewed_uniprot.fasta";
		String diaFilePath = "C:/Users/m334793/Documents/Library/for_context_50perCycle/IT_100ngCurve_100p.dia";
		
		String baseName = diaFilePath.replaceFirst("\\.dia$", "");
		
		// Mass list file 
		String massListPath = "C:/Users/m334793/Documents/Library/for_context_50perCycle/assay.csv";
		
		File fasta = new File(fastaPath);
		File diaFile = new File(diaFilePath);
		File library = new File(libraryPath);
		
//		ArrayList<ScoredFeature> referenceFeatures = new ArrayList<>();
//		ArrayList<ScoredFeature> backgroundFeatures = new ArrayList<>();
	
//		ArrayList<IsolationWindow> targetWindows = IsolationWindowReader.parseMassList(massListPath);
		SearchParameters params = SearchParameterParser.getDefaultParametersObject();
		
		
		// Score features in the .dia file against the library, split the results
		try {
			ContextFeatureScorer.scoreFeatures(library, diaFile, fasta, baseName, massListPath); // run this if the feature file hasn't been processed yet
			File referenceFeatureFile = new File(baseName + "_reference_features.txt");
			File backgroundFeatureFile = new File(baseName + "_background_features.txt");
			
			MProphetExecutionData backgroundData = makeMProphetExecutionData(backgroundFeatureFile, fasta, params, ".background.model.txt");
			MProphetExecutionData referenceData = makeMProphetExecutionData(referenceFeatureFile, fasta, params, ".reference.model.txt");
			
			float peptideFDRThreshold = 0.01f;
			int seed = 1;
			int round = 1;
			
			MProphetResult backgroundMProphetResult = MProphetReiter.executeMProphetTSV(backgroundData, peptideFDRThreshold, seed, params.getAAConstants(), round);
			LinearDiscriminantAnalysis backgroundLDA = backgroundMProphetResult.getLDA();
			
 // 	Use the background LDA model on the reference feature file without retraining
			MProphetResult referenceMProphetResult = MProphetReiter.executeMProphetTSVWithModel(referenceData, peptideFDRThreshold, backgroundLDA, params.getAAConstants());
		
			System.out.println("Finished training lda model on background features. Parsed reference features from " + backgroundFeatureFile.getAbsolutePath());
//			System.out.println("Background passing peptides: " + backgroundMProphetResult.getPassingPeptides().size());
			System.out.println("Finished scoring peptides with background-trained lda model. "
					+ "\nReference passing peptides: " + referenceMProphetResult.getPassingPeptides().size());
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	private static MProphetExecutionData makeMProphetExecutionData(File inputFeatureFile, File fasta, SearchParameters params, String outputSuffix) {

		File peptideOutputFile = new File(inputFeatureFile.getAbsolutePath() + outputSuffix + ".output.txt");
		File peptideDecoyFile = new File(inputFeatureFile.getAbsolutePath() + outputSuffix + ".decoy.txt");

		return new MProphetExecutionData(
				inputFeatureFile,
				fasta,
				peptideOutputFile,
				peptideDecoyFile,
				params
		);
	}

}
