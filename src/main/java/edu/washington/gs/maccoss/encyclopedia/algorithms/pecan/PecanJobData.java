package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public class PecanJobData extends QuantitativeSearchJobData {
	public static final String OUTPUT_FILE_SUFFIX=".pecan.txt";
	public static final String DECOY_FILE_SUFFIX=".pecan.decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".pecan.protein.txt";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".pecan.protein_decoy.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final Optional<ArrayList<FastaPeptideEntry>> targetList;
	private final File fastaFile;
	private final PecanScoringFactory taskFactory;

	public PecanJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, PecanScoringFactory taskFactory) {
		this(targetList, diaFile, fastaFile, getPercolatorExecutionData(diaFile, fastaFile, taskFactory.getParameters()), taskFactory);
	}

	public PecanJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, File outputFile, PecanScoringFactory taskFactory) {
		this(targetList, diaFile, fastaFile, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory);
	}

	public PecanJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, PercolatorExecutionData percolatorFiles, PecanScoringFactory taskFactory) {
		this(targetList, diaFile, null, fastaFile, percolatorFiles, taskFactory);
	}

	public PecanJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, StripeFileInterface diaFileReader, File fastaFile, PercolatorExecutionData percolatorFiles,
			PecanScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, taskFactory.getParameters(), ProgramType.getGlobalVersion().toString());

		this.targetList=targetList;
		this.fastaFile=fastaFile;
		this.taskFactory=taskFactory;
	}

	public static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fastaFile, SearchParameters parameters) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fastaFile,
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}

	private static String getPrefixFromOutput(File outputFile) {
		final String absolutePath = outputFile.getAbsolutePath();

		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - OUTPUT_FILE_SUFFIX.length());
		} else {
			return absolutePath;
		}
	}

	public Optional<ArrayList<FastaPeptideEntry>> getTargetList() {
		return targetList;
	}

	public File getFastaFile() {
		return fastaFile;
	}

	public PecanScoringFactory getTaskFactory() {
		return taskFactory;
	}
	
	@Override
	public String getSearchType() {
		return "Pecan";
	}

	public static String getOutputAbsolutePathPrefix(String absolutePath) {
		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			absolutePath=absolutePath.substring(0, absolutePath.length()-OUTPUT_FILE_SUFFIX.length());
		}
		return absolutePath;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.ELIB);
	}
}