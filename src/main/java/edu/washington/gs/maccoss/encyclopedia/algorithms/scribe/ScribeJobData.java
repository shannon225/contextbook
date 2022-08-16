package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.DDASearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibrarySearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public class ScribeJobData extends DDASearchJobData implements LibrarySearchJobData {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String OUTPUT_FILE_SUFFIX=".scribe.txt";
	public static final String DECOY_FILE_SUFFIX=".scribe.decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".scribe.protein.txt";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".scribe.protein_decoy.txt";
	public static final String FEATURE_FILE_SUFFIX=".scribe.features.txt";

	private final File fastaFile;
	private final LibraryInterface library;
	private final ScribeScoringFactory taskFactory;

	// gui
	public ScribeJobData(File diaFile, File fastaFile, LibraryInterface library, ScribeScoringFactory taskFactory) {
		this(diaFile, fastaFile, library, getPercolatorExecutionData(diaFile, fastaFile, taskFactory.getParameters()), taskFactory);
	}

	// command line
	public ScribeJobData(File diaFile, File fastaFile, LibraryInterface library, File outputFile, ScribeScoringFactory taskFactory) {
		this(diaFile, fastaFile, library, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory);
	}

	// internal only
	private ScribeJobData(File diaFile, File fastaFile, LibraryInterface library, PercolatorExecutionData percolatorFiles, ScribeScoringFactory taskFactory) {
		this(diaFile, null, fastaFile, library, percolatorFiles, taskFactory);
	}

	// used by testing
	public ScribeJobData(File diaFile, StripeFileInterface diaFileReader, File fastaFile, LibraryInterface library, PercolatorExecutionData percolatorFiles,
			ScribeScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, taskFactory.getParameters(), ProgramType.getGlobalVersion().toString());

		this.library=library;
		this.fastaFile=fastaFile;
		this.taskFactory=taskFactory;
	}

	public static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fastaFile, SearchParameters parameters) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fastaFile, 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}

	protected static String getPrefixFromOutput(File outputFile) {
		final String absolutePath = outputFile.getAbsolutePath();

		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - OUTPUT_FILE_SUFFIX.length());
		} else {
			return absolutePath;
		}
	}

	public File getFastaFile() {
		return fastaFile;
	}

	public ScribeScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.DLIB);
	}

	@Override
	public String getSearchType() {
		return "Scribe";
	}
	public LibraryInterface getLibrary() {
		return library;
	}
	
	@Override
	public String getPrimaryScoreName() {
		return taskFactory.getPrimaryScoreName();
	}
}