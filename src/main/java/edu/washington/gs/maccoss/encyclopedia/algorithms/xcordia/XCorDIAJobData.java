package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public class XCorDIAJobData extends QuantitativeSearchJobData {
	public static final String OUTPUT_FILE_SUFFIX=".xcordia.txt";
	public static final String DECOY_FILE_SUFFIX=".xcordia.decoy.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final Optional<ArrayList<FastaPeptideEntry>> targetList;
	private final File fastaFile;
	private final XCorDIAOneScoringFactory taskFactory;

	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, XCorDIAOneScoringFactory taskFactory) {
		this(
				targetList,
				diaFile,
				null,
				fastaFile,
				featureFile,
				outputFile,
				taskFactory
		);
	}

	public XCorDIAJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, StripeFileInterface diaFileReader, File fastaFile, File featureFile, File outputFile, XCorDIAOneScoringFactory taskFactory) {
		super(
				diaFile,
				diaFileReader,
				featureFile,
				outputFile,
				getDecoyFileFromOutput(outputFile),
				taskFactory.getParameters(),
				taskFactory.getVersion()
		);

		this.targetList = targetList;
		this.fastaFile = fastaFile;
		this.taskFactory = taskFactory;
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
			return outputFile.getAbsolutePath();
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
		String absolutePath = getPrefixFromOutput(getOutputFile());
		return new File(absolutePath+LibraryFile.ELIB);
	}

	@Override
	public String getSearchType() {
		return "XCorDIA";
	}
}