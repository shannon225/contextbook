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
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.OutputMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ExternalExecutor;

public class PercolatorExecutor extends ExternalExecutor {
	public static final String PI_0_TAG="pi_0=";
	public static final String V3_01="v3-01";
	public static final String V2_10="v2-10";
	public static final byte DEFAULT_VERSION_NUMBER=3;

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

	PercolatorExecutor(int percolatorVersion, PercolatorExecutionData commandData) {
		super(generateCommand(percolatorVersion, commandData));
	}

	public static Pair<ArrayList<PercolatorPeptide>, Float> executePercolatorTSV(int percolatorVersion, PercolatorExecutionData commandData, float threshold, AminoAcidConstants aaConstants) throws IOException, FileNotFoundException, UnsupportedEncodingException, InterruptedException {
		PercolatorExecutor e=new PercolatorExecutor(percolatorVersion, commandData);
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

	static String[] generateCommand(int percolatorVersion, PercolatorExecutionData commandData) {
		File percolator=getPercolator(percolatorVersion);

		if (percolatorVersion==2) {
			return new String[] {
					percolator.getAbsolutePath(),
					"--results-peptides", commandData.getPeptideOutputFile().getAbsolutePath(),
					"--decoy-results-peptides", commandData.getPeptideDecoyFile().getAbsolutePath(),
					commandData.isUseMinMax()?"-y":"-Y",
					commandData.getInputTSV().getAbsolutePath()
			};
		} else {
			return new String[] {
					percolator.getAbsolutePath(),
					"--results-peptides", commandData.getPeptideOutputFile().getAbsolutePath(),
					"--decoy-results-peptides", commandData.getPeptideDecoyFile().getAbsolutePath(),
					// Commented params below removed when Percolator protein FDR calculations were abandoned
					//"-P", LibraryEntry.DECOY_STRING,
					//"-f", getFastaPlusDecoyFile(commandData.getFastaFile(), commandData.getParameters()).getAbsolutePath(),
					//"--results-proteins", commandData.getProteinOutputFile().getAbsolutePath(),
					//"--decoy-results-proteins", commandData.getProteinDecoyFile().getAbsolutePath(),
					//"--protein-enzyme", commandData.getParameters().getEnzyme().getPercolatorName(),
					//"-g",
					commandData.isUseMinMax()?"-y":"-Y",
					"--no-terminate",
					"-N", Integer.toString(commandData.getParameters().getPercolatorTrainingSetSize()),
					// Note that in EncyclopeDIA 0.9.4 and earlier we did not pass an FDR threshold/test
					// FDR (-t/--testFDR) to Percolator; this meant we always used the default setting (0.01)
					// to evaluate cross-validation results (but the FDR threshold was applied by EncyclopeDIA
					// externally, using a correction for the presence of extra decoys).
					// We now pass this value to Percolator from the search parameters, which simplifies behavior
					// when the training set threshold (-F/--trainFDR) is set to fall back to the test FDR.
					"--testFDR", Float.toString(commandData.getParameters().getPercolatorTestThreshold()),
					"--trainFDR", Float.toString(commandData.getParameters().getPercolatorTrainingSetThreshold()),
					commandData.getInputTSV().getAbsolutePath()
			};
		}
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

	static File getPercolator(int percolatorVersionNumber) {
		String percolatorVersion;
		if (percolatorVersionNumber==2) {
			percolatorVersion=V2_10;
		} else {
			percolatorVersion=V3_01;
		}

		try {
			File percolator=File.createTempFile("Percolator-" + percolatorVersion + "-", ".exe");
			percolator.deleteOnExit();

			OS os=OSDetector.getOS();
			switch (os) {
				case WINDOWS: {
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-"+percolatorVersion+".exe");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);

					// not necessary for the crux version of percolator
					//loadLibraryFile(percolator, "xerces-c_3_1.dll");
					//loadLibraryFile(percolator, "msvcr120.dll");
					//loadLibraryFile(percolator, "msvcp120.dll");

					return percolator;
				}
				case MAC: {
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-"+percolatorVersion+".mac");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					return percolator;
				}
				case LINUX:
					InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/percolator-"+percolatorVersion+".lin");
					Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
					percolator.setExecutable(true);
					return percolator;
			}
			throw new EncyclopediaException("Sorry, Percolator for "+OSDetector.getOSName(os)+" is not set up yet!");
		} catch (IOException ioe) {
			throw new EncyclopediaException("Unexpected exception finding Percolator", ioe);
		}
	}

	static void loadLibraryFile(File percolator, String target) throws IOException {
		File file=new File(percolator.getParentFile(), target);
		file.deleteOnExit();
		InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/"+target);
		Files.copy(is, file.toPath(), StandardCopyOption.REPLACE_EXISTING);
	}
}