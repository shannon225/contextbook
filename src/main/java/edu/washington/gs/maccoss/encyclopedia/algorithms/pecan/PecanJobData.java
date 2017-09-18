package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AbstractSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

public class PecanJobData extends AbstractSearchJobData {
	public static final String OUTPUT_FILE_SUFFIX=".pecan.txt";
	public static final String DECOY_FILE_SUFFIX=".pecan.decoy.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final Optional<ArrayList<FastaPeptideEntry>> targetList;
	private final File fastaFile;
	private final PecanScoringFactory taskFactory;

	public PecanJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, PecanScoringFactory taskFactory) {
		super(diaFile, featureFile, outputFile, new File(getOutputAbsolutePathPrefix(outputFile.getAbsolutePath())+DECOY_FILE_SUFFIX), taskFactory.getParameters(), taskFactory.getVersion());
		this.targetList=targetList;
		this.fastaFile=fastaFile;
		this.taskFactory=taskFactory;
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
}