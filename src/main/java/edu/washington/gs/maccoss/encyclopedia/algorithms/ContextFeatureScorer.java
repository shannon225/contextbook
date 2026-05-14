package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.EncyclopediaTwo;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaTwoJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.context.IsolationWindow;
import edu.washington.gs.maccoss.encyclopedia.context.IsolationWindowReader;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class ContextFeatureScorer {

	public static void main(String[] args) throws IOException, SQLException, InterruptedException, DataFormatException {
		// Inputs for Search
		String rawFilePath = "C:/Users/m334793/Documents/Library/cd14_combined.dia"; // Raw file to search
		String libraryFilePath = "C:/Users/m334793/Documents/Library/easyspray_lit_immune_library.elib"; // Library to
																											// search
																											// against
		String fastaPath = "C:/Users/m334793/Documents/Library/human_uniprot_2025dec12.fasta"; // fasta file for serach
		String baseName = rawFilePath.replaceFirst("\\.dia$", "");

		// Mass list file
		String massListPath = "C:/Users/m334793/Documents/Library/assay7.csv";

		final File fasta = new File(fastaPath);
		File rawFile = new File(rawFilePath);
		File library = new File(libraryFilePath);

		try {
			ArrayList<ScoredFeature> bestFeature = scoreFeatures(library, rawFile, fasta, baseName, massListPath);
			System.out.println(bestFeature.get(0).getSequence());
			System.out.println(bestFeature.get(0).getPrimary());

		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	private static boolean isFeatureOnMassList(
	        ScoredFeature feature,
	        ArrayList<IsolationWindow> targetWindows,
	        double halfWindowWidth) {

	    double featureMz = feature.getMz();

	    for (IsolationWindow window : targetWindows) {
	        double targetMz = window.getTargetMz();
	        double mzStart = targetMz - halfWindowWidth;
	        double mzStop = targetMz + halfWindowWidth;

	        if (featureMz >= mzStart && featureMz <= mzStop) {
	            return true;
	        }
	    }

	    return false;
	}

	private static void writeScoredFeatures(File outputFile, ArrayList<ScoredFeature> features, String header)
	        throws IOException {
	    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
	        writer.write(header);
	        writer.newLine();

	        for (ScoredFeature feature : features) {
	            writer.write(feature.getOriginalLine());
	            writer.newLine();
	        }
	    }
	}

	public static ArrayList<ScoredFeature> scoreFeatures(File library, File rawFile, File fasta, String baseName,
			String massListPath) throws IOException, SQLException, DataFormatException, InterruptedException {

		// Run an Encyclopedia job
		SearchParameters params = SearchParameterParser.getDefaultParametersObject();
		LibraryScoringFactory scoringForLibrary = EncyclopediaScoringFactory.getDefaultScoringFactory(params);
		LibraryInterface interfaceForLibrary = BlibToLibraryConverter.getFile(library, fasta, params);

		EncyclopediaTwoJobData job = new EncyclopediaTwoJobData(rawFile, fasta, interfaceForLibrary,
				interfaceForLibrary, rawFile, scoringForLibrary);

		ProgressIndicator progress = new EmptyProgressIndicator(true);
		StripeFileInterface interfaceForStripeFile = job.getDiaFileReader();

		// Run Encyclopedia job to get the feature file
		File featuresToSplit = job.getPercolatorFiles().getInputTSV();

		if (featuresToSplit.exists() && featuresToSplit.canRead()) {
			System.out.println("Feature file already exists, skipping feature calculation!");
			System.out.println(featuresToSplit.getAbsolutePath());
		} else {
			System.out.println("Calculating features...");
			EncyclopediaTwo.generateFeatureFile(progress, interfaceForLibrary, job,
					interfaceForStripeFile, java.util.Optional.empty());
		}

//		String outputPathForUniqueFeatures = baseName + "all_features.txt";
		ArrayList<ScoredFeature> uniqueFeatures = new ArrayList<>();
		ArrayList<ScoredFeature> uniqueFeaturesList = uniqueFeatures;
		HashMap<String, ScoredFeature> bestFeatureByPeptide = new HashMap<>();
		String header;

		// read all rows the feature file
		try (BufferedReader br = new BufferedReader(new FileReader(featuresToSplit))) {
			header = br.readLine();

			String line;
			while ((line = br.readLine()) != null) {
				String columns[] = line.split("\t", -1);

				double mz = Double.parseDouble(columns[27]);
				boolean isDecoy = Integer.parseInt(columns[1]) == -1;
				float primary = Float.parseFloat(columns[3]);
				float retentionTime = Float.parseFloat(columns[29]);
				String sequence = columns[30];
				String protein = columns[31];

				ScoredFeature feature = new ScoredFeature(mz, isDecoy, primary, retentionTime, sequence, protein, line);

				uniqueFeaturesList.add(feature);

				ScoredFeature currentBest = bestFeatureByPeptide.get(sequence);

				// Take the peptide with a higher primary score and place it on a new list

				if (currentBest == null || feature.getPrimary() > currentBest.getPrimary()) {
					bestFeatureByPeptide.put(sequence, feature);
				}

			}
		}
		ArrayList<ScoredFeature> bestFeatures = new ArrayList<>(bestFeatureByPeptide.values());
		bestFeatures.sort(Comparator.comparing(ScoredFeature::getPrimary).reversed());

		System.out.println("Unique scored features have been found " + bestFeatures.size());


		// Output Paths
		String referenceOutputPath = baseName + "_reference_features.txt";
		String backgroundOutputPath = baseName + "_background_features.txt";
		String referenceDecoyOutputPath = baseName + "_reference_decoy_features.txt";
		String backgroundDecoyOutputPath = baseName + "_background_decoy_features.txt";

		// Output Files
		File referenceOutput = new File(referenceOutputPath);
		File backgroundOutput = new File(backgroundOutputPath);
		File referenceDecoyOutput = new File(referenceDecoyOutputPath);
		File backgroundDecoyOutput = new File(backgroundDecoyOutputPath);

		// Target mass list
		ArrayList<IsolationWindow> targetWindows = IsolationWindowReader.parseMassList(massListPath);
		System.out.println(targetWindows.size() + " windows cataloged from the mass list");

		ArrayList<ScoredFeature> referenceFeatures = new ArrayList<>();
		ArrayList<ScoredFeature> backgroundFeatures = new ArrayList<>();
		ArrayList<ScoredFeature> referenceDecoyFeatures = new ArrayList<>();
		ArrayList<ScoredFeature> backgroundDecoyFeatures = new ArrayList<>();
		
		double halfWindowWidth = 0.35;

		for (ScoredFeature feature : bestFeatures) {
		    boolean isOnMassList = isFeatureOnMassList(feature, targetWindows, halfWindowWidth);

		    if (!feature.isDecoy() && isOnMassList) {
		        referenceFeatures.add(feature);
		    } else if (!feature.isDecoy() && !isOnMassList) {
		        backgroundFeatures.add(feature);
		    } else if (feature.isDecoy() && isOnMassList) {
		        referenceDecoyFeatures.add(feature);
		    } else {
		        backgroundDecoyFeatures.add(feature);
		    }
		}
		writeScoredFeatures(referenceOutput, referenceFeatures, header);
		writeScoredFeatures(backgroundOutput, backgroundFeatures, header);
		writeScoredFeatures(referenceDecoyOutput, referenceDecoyFeatures, header);
		writeScoredFeatures(backgroundDecoyOutput, backgroundDecoyFeatures, header);

		System.out.println("Reference target features: " + referenceFeatures.size());
		System.out.println("Background target features: " + backgroundFeatures.size());
		System.out.println("Reference decoy features: " + referenceDecoyFeatures.size());
		System.out.println("Background decoy features: " + backgroundDecoyFeatures.size());
	return bestFeatures;
}

	public static void partitionFeatureFile() {

	}
}
