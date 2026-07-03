package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class MProphetBatchRunner {

	public static final Pattern MASKED_PATTERN = Pattern.compile("masked(\\d+)");
	
	public static void mainRunOneFileMultipleTimes(String[] args) throws Exception {
		if (args.length != 4) {
			Logger.errorLine("MProphetReiterMultiSeed requires four parameters in order:");
			Logger.logLine("  1) Input TSV");
			Logger.logLine("  2) Input FASTA");
			Logger.logLine("  3) Threshold (e.g., 0.01)");
			Logger.logLine("  4) Number of seeds (e.g., 10 to run with seeds 1...10)");
			System.exit(1);
		}

		// Parse CLI arguments (same inputs as existing code, but 4th = number of seeds)
		File inputTSV = new File(args[0]);
		File fastaFile = new File(args[1]);
		float threshold = Float.parseFloat(args[2]);
		int numSeeds = Integer.parseInt(args[3]);

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();

		// Determine base directory and create designated output folder
		File baseDir = inputTSV.getParentFile();
		if (baseDir == null) {
			baseDir = new File(".");
		}

		File outputDir = new File(baseDir, "MProphet_Encyclopedia_Output");
		if (!outputDir.exists() && !outputDir.mkdirs()) {
			throw new IOException("Could not create output directory: " + outputDir.getAbsolutePath());
		}

		String baseName = stripExtension(inputTSV.getName());

		try {
			Logger.logLine("Running MProphetReiter in multi-seed mode.");
			Logger.logLine("Input TSV: " + inputTSV.getAbsolutePath());
			Logger.logLine("Input FASTA: " + fastaFile.getAbsolutePath());
			Logger.logLine("Peptide FDR threshold: " + threshold);
			Logger.logLine("Number of seeds: " + numSeeds);
			Logger.logLine("Output directory: " + outputDir.getAbsolutePath());

			for (int seed = 1; seed <= numSeeds; seed++) {
				// Per-seed output files, in the designated folder
				File peptideOutputFile = new File(outputDir, baseName + "_seed" + seed + ".output.txt");
				File peptideDecoyFile = new File(outputDir, baseName + "_seed" + seed + ".decoy.txt");

				Logger.logLine("");
				Logger.logLine("Running MProphet for seed " + seed);
				Logger.logLine("Target output: " + peptideOutputFile.getAbsolutePath());
				Logger.logLine("Decoy output:  " + peptideDecoyFile.getAbsolutePath());

				MProphetExecutionData execData = new MProphetExecutionData(inputTSV, fastaFile, peptideOutputFile,
						peptideDecoyFile, params);

				MProphetReiter.executeMProphetTSV(execData, threshold, seed, params.getAAConstants(), 1);
			}

			Logger.logLine("Multi-seed MProphetReiter run completed successfully.");
		} catch (Exception e) {
			e.printStackTrace();
			throw e;

		} finally {
			
		}
	}

	public static void main(String[] args) throws Exception {
		if (args.length != 4) {
			Logger.errorLine("MProphetReiterMultiSeed requires four parameters in order:");
			Logger.logLine("  1) Input Folder");
			Logger.logLine("  2) Input FASTA");
			Logger.logLine("  3) Threshold (e.g., 0.01)");
	//		Logger.logLine("  4) Number of seeds (e.g., 10 to run with seeds 1...10)");
			System.exit(1);
		}

		// Parse CLI arguments (same inputs as existing code, but 4th = number of seeds)
		File inputDir = new File(args[0]);
		File fastaFile = new File(args[1]);
		float threshold = Float.parseFloat(args[2]);
//		int numSeeds = Integer.parseInt(args[3]);

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();

		// Determine base directory and create designated output folder
		if (!inputDir.exists() || !inputDir.isDirectory()) {
		    throw new IOException("Input location is not a directory: " + inputDir.getAbsolutePath());
		}
		
//		File baseDir = inputDir.getParentFile();
//		if (baseDir == null) {
//			baseDir = new File(".");
//		}

		File outputDir = new File(inputDir, "mProphet_output");
		if (!outputDir.exists() && !outputDir.mkdirs()) {
			throw new IOException("Could not create output directory: " + outputDir.getAbsolutePath());
		}

		File[] inputFiles = inputDir.listFiles((dir, name) -> name.endsWith(".tsv") || name.endsWith(".features.txt"));
		Arrays.sort(inputFiles, Comparator.comparing(File::getName));
		
		try {
			Logger.logLine("Running MProphetReiter in multi-seed mode.");
			Logger.logLine("Input Folder: " + inputDir.getAbsolutePath());
			Logger.logLine("Input FASTA: " + fastaFile.getAbsolutePath());
			Logger.logLine("Peptide FDR threshold: " + threshold);
//			Logger.logLine("Number of seeds: " + numSeeds);
			Logger.logLine("Output directory: " + outputDir.getAbsolutePath());

			for (File inputTSV : inputFiles) {
				int seed = extractSeedFromFileName(inputTSV);
				String baseName = stripExtension(inputTSV.getName());

				// Per-seed output files, in the designated folder
				File peptideOutputFile = new File(outputDir, baseName + ".output.txt");
				File peptideDecoyFile = new File(outputDir, baseName + ".decoy.txt");

				Logger.logLine("");
				Logger.logLine("Running MProphet for file " + baseName);
				Logger.logLine("Target output: " + peptideOutputFile.getAbsolutePath());
				Logger.logLine("Decoy output:  " + peptideDecoyFile.getAbsolutePath());

				MProphetExecutionData execData = new MProphetExecutionData(inputTSV, fastaFile, peptideOutputFile,
						peptideDecoyFile, params);

				MProphetReiter.executeMProphetTSV(execData, threshold, seed, params.getAAConstants(), 1);
			}

			Logger.logLine("Multi-seed MProphetReiter run completed successfully.");
		} catch (Exception e) {
			e.printStackTrace();
			throw e;

		} finally {
			
		}
	}
	private static String stripExtension(String name) {
		int dot = name.lastIndexOf('.');
		return (dot == -1) ? name : name.substring(0, dot);
	}

	private static int extractSeedFromFileName(File file) {
	    Matcher matcher = MASKED_PATTERN.matcher(file.getName());

	    if (!matcher.find()) {
	        throw new IllegalArgumentException(
	            "Could not extract maskedN seed from filename: " + file.getName()
	        );
	    }

	    return Integer.parseInt(matcher.group(1));
	}
}
