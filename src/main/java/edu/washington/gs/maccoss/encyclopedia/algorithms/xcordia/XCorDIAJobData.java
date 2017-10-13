package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public class XCorDIAJobData extends QuantitativeSearchJobData {
	public static final String OUTPUT_FILE_SUFFIX=".xcordia.txt";
	public static final String DECOY_FILE_SUFFIX=".xcordia.decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".xcordia.protein.txt";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".xcordia.protein_decoy.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final Optional<ArrayList<FastaPeptideEntry>> targetList;
	private final File fastaFile;
	private final XCorDIAOneScoringFactory taskFactory;

	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, XCorDIAOneScoringFactory taskFactory) {
		this(targetList, diaFile, fastaFile, getPercolatorExecutionData(diaFile, fastaFile), taskFactory);
	}

	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, File outputFile, XCorDIAOneScoringFactory taskFactory) {
		this(targetList, diaFile, fastaFile, getPercolatorExecutionData(outputFile, fastaFile), taskFactory);
	}

	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, PercolatorExecutionData percolatorFiles, XCorDIAOneScoringFactory taskFactory) {
		this(targetList, diaFile, null, fastaFile, percolatorFiles, taskFactory);
	}

	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, StripeFileInterface diaFileReader, File fastaFile, PercolatorExecutionData percolatorFiles,
			XCorDIAOneScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, taskFactory.getParameters(), taskFactory.getVersion());

		this.targetList=targetList;
		this.fastaFile=fastaFile;
		this.taskFactory=taskFactory;
	}

	private static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fastaFile) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fastaFile, 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX));
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

	public XCorDIAOneScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.ELIB);
	}

	@Override
	public String getSearchType() {
		return "XCorDIA";
	}
}