package edu.washington.gs.maccoss.encyclopedia.context;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetReiter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;

/**
 * Context-mode mProphet driver.
 *
 * End-to-end: scores a single PRM file against a library, splits features
 * into reference (peptides in the mass list) vs background (everything else),
 * trains an LDA mProphet model on the background features, and applies that
 * model to the reference features.
 *
 * Outputs (written alongside each split feature file):
 *   <basename>_reference.features.pep.output.txt   (target peptides)
 *   <basename>_reference.features.pep.decoy.txt    (decoy peptides)
 *   <basename>_background.features.pep.output.txt
 *   <basename>_background.features.pep.decoy.txt
 * Plus diagnostic plots in -plotsdir.
 */
public class ContextMProphetExecutor {

	private static final float  DEFAULT_FDR  = 0.01f;
	private static final int    DEFAULT_SEED = 1;
	private static final String DEFAULT_PLOTS_DIR = "mprophet_plots";

	public static void main(String[] args) throws IOException {
		HashMap<String, String> arguments = CommandLineParser.parseArguments(args);

		if (args.length == 0
				|| arguments.containsKey("-h")
				|| arguments.containsKey("-help")
				|| arguments.containsKey("--help")) {
			printHelp();
			System.exit(args.length == 0 ? 1 : 0);
		}

		File dia       = requiredFile(arguments, "-i");
		File library   = requiredFile(arguments, "-l");
		File fasta     = requiredFile(arguments, "-f");
		File massList  = requiredFile(arguments, "-massList");
		float fdr      = Float.parseFloat(arguments.getOrDefault("-fdr",  Float.toString(DEFAULT_FDR)));
		int   seed     = Integer.parseInt(arguments.getOrDefault("-seed", Integer.toString(DEFAULT_SEED)));
		File  plotsDir = new File(arguments.getOrDefault("-plotsdir", DEFAULT_PLOTS_DIR));

		try {
			runEndToEnd(library, fasta, dia, massList, fdr, seed, plotsDir);
		} catch (Exception e) {
			Logger.errorLine("ContextMProphetExecutor failed: " + e.getMessage());
			e.printStackTrace();
			System.exit(2);
		}
	}

	/**
	 * Full pipeline for one PRM file: score -> split -> train(bg) -> apply(ref).
	 */
	public static void runEndToEnd(File library, File fasta, File dia, File massList,
	                               float fdr, int seed, File plotsDir) throws Exception {

		if (!dia.exists())      throw new IOException("Input file not found: " + dia);
		if (!library.exists())  throw new IOException("Library file not found: " + library);
		if (!fasta.exists())    throw new IOException("FASTA file not found: " + fasta);
		if (!massList.exists()) throw new IOException("Mass list file not found: " + massList);

		// Strip .dia or .mzML (case-insensitive). EncyclopeDIA accepts both.
		String baseName = dia.getAbsolutePath().replaceFirst("(?i)\\.(dia|mzML)$", "");

		Logger.logLine("ContextFeatureScorer: scoring " + dia.getName());
		ContextFeatureScorer.scoreFeatures(library, dia, fasta, baseName, massList.getAbsolutePath());

		File referenceFeatures  = new File(baseName + "_reference.features.txt");
		File backgroundFeatures = new File(baseName + "_background.features.txt");

		runTrainApply(referenceFeatures, backgroundFeatures, fasta, fdr, seed, plotsDir);
	}

	/**
	 * Train mProphet LDA on {@code backgroundFeatures}, apply to {@code referenceFeatures}.
	 * Exposed publicly so future rescorers can reuse the same split feature files.
	 */
	public static void runTrainApply(File referenceFeatures, File backgroundFeatures,
	                                 File fasta, float fdr, int seed, File plotsDir) throws Exception {

		SearchParameters params = SearchParameterParser.getDefaultParametersObject();

		MProphetExecutionData backgroundData = buildExecutionData(backgroundFeatures, fasta, params);
		MProphetExecutionData referenceData  = buildExecutionData(referenceFeatures,  fasta, params);

		final int round = 1;

		Logger.logLine("Training mProphet LDA on background features: " + backgroundFeatures.getName());
		MProphetResult backgroundResult = MProphetReiter.executeMProphetTSV(
				backgroundData, fdr, seed, params.getAAConstants(), round);
		LinearDiscriminantAnalysis backgroundLDA = backgroundResult.getLDA();

		Logger.logLine("Applying background-trained LDA to reference features: " + referenceFeatures.getName());
		MProphetResult referenceResult = MProphetReiter.executeMProphetTSVWithModel(
				referenceData, fdr, backgroundLDA, params.getAAConstants());

		Logger.logLine("Reference passing peptides: " + referenceResult.getPassingPeptides().size());

		if (!plotsDir.exists() && !plotsDir.mkdirs()) {
			Logger.errorLine("Could not create plots directory: " + plotsDir);
		}
		ContextMProphetPlotter.plotContextMProphetResults(
				backgroundData.getPeptideOutputFile(),
				backgroundData.getPeptideDecoyFile(),
				referenceData.getPeptideOutputFile(),
				referenceData.getPeptideDecoyFile(),
				plotsDir);
	}

	private static MProphetExecutionData buildExecutionData(File inputFeatures, File fasta, SearchParameters params) {
		String base = inputFeatures.getAbsolutePath().replaceAll("\\.txt$", "");
		File peptideOutputFile = new File(base + ".pep.output.txt");
		File peptideDecoyFile  = new File(base + ".pep.decoy.txt");
		return new MProphetExecutionData(inputFeatures, fasta, peptideOutputFile, peptideDecoyFile, params);
	}

	private static File requiredFile(HashMap<String, String> args, String flag) {
		String value = args.get(flag);
		if (value == null) {
			Logger.errorLine("Missing required argument: " + flag);
			printHelp();
			System.exit(1);
		}
		return new File(value);
	}

	private static void printHelp() {
		Logger.timelessLogLine("ContextMProphetExecutor");
		Logger.timelessLogLine("Context-mode mProphet: score one PRM file, split by mass-list membership,");
		Logger.timelessLogLine("train LDA on background features, apply to reference features.");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Required:");
		Logger.timelessLogLine("  -i        <file>   input .dia (or .mzML)");
		Logger.timelessLogLine("  -l        <file>   library (.elib preferred, .dlib accepted)");
		Logger.timelessLogLine("  -f        <file>   FASTA protein database");
		Logger.timelessLogLine("  -massList <file>   assay / mass-list .txt");
		Logger.timelessLogLine("Optional:");
		Logger.timelessLogLine("  -fdr      <float>  peptide FDR threshold (default: " + DEFAULT_FDR + ")");
		Logger.timelessLogLine("  -seed     <int>    random seed for LDA training (default: " + DEFAULT_SEED + ")");
		Logger.timelessLogLine("  -plotsdir <dir>    diagnostic plot directory (default: " + DEFAULT_PLOTS_DIR + ")");
	}
}
