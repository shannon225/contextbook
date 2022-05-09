package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.filewriters.FastaWriter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.FileConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ExternalExecutor;

public class PercolatorExecutor extends ExternalExecutor {
	public static final String PI_0_TAG="pi_0=";

	/**
	 * Instead use {@link PercolatorVersion#DEFAULT_VERSION}.
	 */
	@Deprecated
	public static final PercolatorVersion DEFAULT_VERSION_NUMBER=PercolatorVersion.DEFAULT_VERSION;

	/**
	 * The default value that Percolator uses for the {@code -t/--testFDR} parameter.
	 *
	 * @see SearchParameters#getPercolatorTestThreshold() for more information about
	 *      this value's use in EncyclopeDIA 0.9.4 and earlier.
	 */
	public static final float PERCOLATOR_DEFAULT_TEST_THRESHOLD = 0.01f;

	/**
	 * The value that Percolator uses for the {@code -F/--trainFDR}
	 * parameter to indicate that the training set FDR should match
	 * the test FDR/peptide FDR threshold (the default behavior).
	 */
	public static final float PERCOLATOR_TRAINING_THRESHOLD_FALLBACK_VALUE = 0f;

	/**
	 * The default value that EncyclopeDIA uses for Percolator's
	 * {@code -N} parameter, which sets the number of PSMs to
	 * use as the training set.
	 */
	public static final int DEFAULT_TRAINING_SET_SIZE = 500000;

	/**
	 * The default value that EncyclopeDIA uses for Percolator's
	 * {@code -F/--trainFDR}, which sets the FDR used to select
	 * positive examples from the training set.
	 * By default, use the value {@link #PERCOLATOR_TRAINING_THRESHOLD_FALLBACK_VALUE}
	 * to use Percolator's default, which is to use the same FDR
	 * as the peptide FDR threshold ({@code -t/--testFDR}). Note
	 * that this fallback may behave in an unexpected way due to
	 * the way that EncyclopeDIA handles peptide FDR thresholds;
	 * see {@link SearchParameters#getPercolatorTestThreshold()}
	 * for more information about the {@code -t/--testFDR} parameter.
	 */
	public static final float DEFAULT_TRAINING_THRESHOLD = PERCOLATOR_TRAINING_THRESHOLD_FALLBACK_VALUE;

	private static final Pattern PERCOLATOR_VERSION_PATTERN = Pattern.compile("Percolator version (.+),");
	private static final String SELECTING_PI_0 = "Selecting pi_0=";
	private static final String ERROR_PREFIX = "Error : ";
	private static final String BAD_ALLOCATION = "bad allocation";
	private static final String EXCEPTION_CAUGHT_PREFIX = "Exception caught: ";

	PercolatorExecutor(PercolatorVersion percolatorVersion, PercolatorExecutionData commandData, int round) {
		super(generateCommand(percolatorVersion, commandData, round));
	}

	public static Pair<ArrayList<PercolatorPeptide>, Float> executePercolatorTSV(PercolatorVersion percolatorVersion, PercolatorExecutionData commandData, float threshold, AminoAcidConstants aaConstants, int round) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		PercolatorExecutor e=new PercolatorExecutor(percolatorVersion, commandData, round);
		BlockingQueue<OutputMessage> result=e.start();

		Float pi0=null;
		String errorMessage=null;
		Optional<String> percolatorExecutableVersion = Optional.empty();
		while (!e.isFinished()||!result.isEmpty()) {
			if (!result.isEmpty()) {
				OutputMessage data=result.take();
				if (!data.isStdOutput()) {
					if (!percolatorExecutableVersion.isPresent()) {
						String message = data.getMessage();
						percolatorExecutableVersion = getPercolatorVersionFromOutput(message);
					}

					Logger.logLine(data.getMessage());
					errorMessage = getErrorMessage(data);

					if (null == errorMessage) {
						final String trim = data.getMessage().trim();
						if (trim.startsWith(SELECTING_PI_0)) {
							try {
								pi0 = Float.parseFloat(trim.substring(SELECTING_PI_0.length()));
							} catch (NumberFormatException nfe) {
								Logger.errorLine("Error parsing pi0 from [" + trim + "]");
							}
						}
					}
				}
			} else {
				Thread.sleep(10);
			}
		}

		if (errorMessage!=null) {
			throw new EncyclopediaException(errorMessage);
		}

		checkResult(e);

		try {
		    Files.write(commandData.getPeptideOutputFile().toPath(), (PI_0_TAG+pi0+System.lineSeparator()).getBytes(), StandardOpenOption.APPEND);
		    // if round 1, then start the model file over, otherwise append weights to model file
		    if (round==1) {
		    		commandData.getWeightsFile(round).renameTo(commandData.getModelFile());
		    } else {
		    		FileConcatenator.saveConcatenatedFile(commandData.getModelFile(), commandData.getWeightsFile(round));
		    		commandData.getWeightsFile(round).delete();
		    }
		}catch (IOException ioe) {
			throw new EncyclopediaException("Error appending to Percolator text file", ioe);
		}

		commandData.setPercolatorExecutableVersion(percolatorExecutableVersion.orElse(null));

		Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides=PercolatorReader.getPassingPeptidesFromTSV(commandData.getPeptideOutputFile(), threshold, aaConstants, false);

		return passingPeptides;
	}

	static String getErrorMessage(OutputMessage data) {
		final String trim=data.getMessage().trim();

		final String errorMessage;
		if (trim.startsWith(ERROR_PREFIX)) {
			errorMessage = trim.substring(ERROR_PREFIX.length());
		} else if (trim.startsWith(EXCEPTION_CAUGHT_PREFIX)) {
			errorMessage = trim.substring(EXCEPTION_CAUGHT_PREFIX.length());
		} else if (trim.contains(BAD_ALLOCATION)) {
			errorMessage = trim;
		} else {
			errorMessage = null;
		}
		return errorMessage;
	}

	static Optional<String> getPercolatorVersionFromOutput(String standardOutputLine) {
		Matcher matcher = PERCOLATOR_VERSION_PATTERN.matcher(standardOutputLine);
		if (matcher.find()) {
			return Optional.ofNullable(matcher.group(1));
		} else {
			return Optional.empty();
		}
	}

	private static void checkResult(PercolatorExecutor e) throws EncyclopediaException {
		if (0 != e.getResultCode()) {
			throw new EncyclopediaException("Percolator exited with non-zero status: " + e.getResultCode());
		}
	}

	static String parsePeptideSequence(String peptideString) {
		return peptideString.substring(peptideString.indexOf('.')+1, peptideString.lastIndexOf('.'));
	}

	static String[] generateCommand(PercolatorVersion percolatorVersion, PercolatorExecutionData commandData, int round) {
		File percolator = percolatorVersion.getPercolator();

		ArrayList<String> params=new ArrayList<>();
		
		params.add(percolator.getAbsolutePath());
		params.add("--results-peptides"); params.add(commandData.getPeptideOutputFile().getAbsolutePath());
		params.add("--weights"); params.add(commandData.getWeightsFile(round).getAbsolutePath());
		params.add("--decoy-results-peptides"); params.add(commandData.getPeptideDecoyFile().getAbsolutePath());
		if (commandData.isUseMinMax()) {
			params.add("-y");
		} else {
			params.add("-Y");
		}
		if (commandData.getPercolatorModelFile().isPresent()&&commandData.getPercolatorModelFile().get().canRead()) {
			File modelFile = commandData.getPercolatorModelFile().get();
			try {
				int actualRound=Math.min(round, FileConcatenator.getNumberOfSubFiles(modelFile));
				File model=FileConcatenator.extractFile(modelFile, actualRound);
				Logger.logLine("Extracting weights from "+modelFile.getName()+" ("+round+","+actualRound+")");
				if (round!=actualRound) {
					Logger.errorLine("Couldn't extract specific model for round "+round+", using last model available (round "+actualRound+")");
				}
				params.add("--init-weights"); params.add(model.getAbsolutePath());
				params.add("--maxiter"); params.add("0");
			} catch (IOException ioe) {
				Logger.errorLine("Problem extracting Percolator weights from "+modelFile.getName()+". Continuing without using weights...");
				Logger.errorException(ioe);
			}
		}
		
		if (percolatorVersion.getMajorVersion()>2) {
			params.add("--no-terminate");
			params.add("-N"); params.add(Integer.toString(commandData.getParameters().getPercolatorTrainingSetSize()));
			params.add("--testFDR"); params.add(Float.toString(commandData.getParameters().getPercolatorTestThreshold()));
			params.add("--trainFDR"); params.add(Float.toString(commandData.getParameters().getPercolatorTrainingSetThreshold()));
		}
		params.add(commandData.getInputTSV().getAbsolutePath());
		
		return params.toArray(new String[params.size()]);
	}

	public static File getFastaPlusDecoyFile(File fasta, SearchParameters parameters) {
		File fastaPlusDecoy=new File(fasta.getParentFile(), parameters.getEnzyme().getPercolatorName()+"_"+fasta.getName());
		if (fastaPlusDecoy.exists()&&fastaPlusDecoy.canRead()) return fastaPlusDecoy;

		Logger.logLine("Generating reverse-concatenated FASTA: "+fastaPlusDecoy.getName());
		FastaWriter writer=new FastaWriter(fastaPlusDecoy);
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(fasta, parameters);
		for (FastaEntryInterface entry : entries) {
			writer.write(entry);
			FastaEntry reverse=new FastaEntry(entry.getFilename(), LibraryEntry.DECOY_STRING+entry.getAnnotation(), parameters.getEnzyme().reverseProtein(entry.getSequence(), parameters.getAAConstants()));
			writer.write(reverse);
		}
		writer.close();

		return fastaPlusDecoy;
	}
}