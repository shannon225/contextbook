package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public class EncyclopediaJobData extends QuantitativeSearchJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".encyclopedia.protein_decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".encyclopedia.protein.txt";
	public static final String DECOY_FILE_SUFFIX=".encyclopedia.decoy.txt";
	public static final String OUTPUT_FILE_SUFFIX=".encyclopedia.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final LibraryInterface library;
	private final LibraryScoringFactory taskFactory;

	public EncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(diaFile, null, getPercolatorExecutionData(diaFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, File outputFile, LibraryScoringFactory taskFactory) {
		this(diaFile, null, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(diaFile, null, percolatorFiles, parameters, version, library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, StripeFileInterface diaFileReader, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, parameters, version);

		this.library = library;
		this.taskFactory = taskFactory;
	}

	private static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fastaFile, SearchParameters parameters) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fastaFile,
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
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
		return new EncyclopediaJobData(getDiaFile(), diaFileReader, getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}

	public LibraryInterface getLibrary() {
		return library;
	}

	public LibraryScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.ELIB);
	}

	@Override
	public String getSearchType() {
		return "EncyclopeDIA";
	}
}