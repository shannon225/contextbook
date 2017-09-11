package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

import java.io.File;

public class EncyclopediaJobData extends QuantitativeSearchJobData {
	public static final String LOG_FILE_SUFFIX = ".log";
	public static final String DECOY_FILE_SUFFIX = ".encyclopedia.decoy.txt";
	public static final String OUTPUT_FILE_SUFFIX = ".encyclopedia.txt";
	public static final String FEATURE_FILE_SUFFIX = ".features.txt";

	private final LibraryInterface library;
	private final LibraryScoringFactory taskFactory;

	public EncyclopediaJobData(File diaFile, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(
				diaFile,
				null,
				getFeatureFile(diaFile),
				getOutputFile(diaFile),
				taskFactory.getParameters(),
				taskFactory.getVersion(),
				library,
				taskFactory
		);
	}

	public EncyclopediaJobData(File diaFile, LibraryInterface library, File outputFile, LibraryScoringFactory taskFactory) {
		this(
				diaFile,
				null,
				getFeatureFileFromOutput(outputFile),
				outputFile,
				taskFactory.getParameters(),
				taskFactory.getVersion(),
				library,
				taskFactory
		);
	}

	public EncyclopediaJobData(File diaFile, File featureFile, File outputFile, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(
				diaFile,
				null,
				featureFile,
				outputFile,
				parameters,
				version,
				library,
				taskFactory
		);
	}

	public EncyclopediaJobData(File diaFile, StripeFileInterface diaFileReader, LibraryInterface library, File outputFile, LibraryScoringFactory taskFactory) {
		this(
				diaFile,
				diaFileReader,
				getFeatureFileFromOutput(outputFile),
				outputFile,
				taskFactory.getParameters(),
				taskFactory.getVersion(),
				library,
				taskFactory
		);
	}

	public EncyclopediaJobData(File diaFile, StripeFileInterface diaFileReader, File featureFile, File outputFile, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, diaFileReader, featureFile, outputFile, getDecoyFileFromOutput(outputFile), parameters, version);

		this.library = library;
		this.taskFactory = taskFactory;
	}

	static File getFeatureFile(File diaFile) {
		return new File(diaFile.getAbsolutePath() + FEATURE_FILE_SUFFIX);
	}

	static File getFeatureFileFromOutput(File outputFile) {
		return getFeatureFile(new File(getPrefixFromOutput(outputFile)));
	}

	static File getOutputFile(File diaFile) {
		return new File(diaFile.getAbsolutePath() + OUTPUT_FILE_SUFFIX);
	}

	static File getDecoyFile(File diaFile) {
		return new File(diaFile.getAbsolutePath() + DECOY_FILE_SUFFIX);
	}

	static File getDecoyFileFromOutput(File outputFile) {
		return getDecoyFile(new File(getPrefixFromOutput(outputFile)));
	}

	static String getPrefixFromOutput(File outputFile) {
		final String absolutePath = outputFile.getAbsolutePath();

		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - OUTPUT_FILE_SUFFIX.length());
		} else {
			return absolutePath;
		}
	}

	public EncyclopediaJobData updateTaskFactory(LibraryScoringFactory taskFactory) {
		return new EncyclopediaJobData(getDiaFile(), diaFileReader, getFeatureFile(), getOutputFile(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}

	public LibraryInterface getLibrary() {
		return library;
	}

	public LibraryScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getFirstPassPercolator() {
		String absolutePath = getPrefixFromOutput(getOutputFile());
		return new File(absolutePath + ".first_round.txt");
	}

	public File getFirstPassPercolatorDecoy() {
		String absolutePath = getPrefixFromOutput(getOutputFile());
		return new File(absolutePath + ".first_round.decoy.txt");
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getOutputFile());
		return new File(absolutePath + LibraryFile.ELIB);
	}

	@Override
	public String getSearchType() {
		return "EncyclopeDIA";
	}
}