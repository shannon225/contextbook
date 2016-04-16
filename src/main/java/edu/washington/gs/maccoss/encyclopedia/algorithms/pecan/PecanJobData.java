package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;

public class PecanJobData extends SearchJobData {
	private final Optional<ArrayList<FastaPeptideEntry>> targetList;
	private final File fastaFile;
	private final PecanScoringFactory taskFactory;

	public PecanJobData(Optional<ArrayList<FastaPeptideEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, PecanScoringFactory taskFactory) {
		super(diaFile, featureFile, outputFile, taskFactory.getParameters(), taskFactory.getVersion());
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
}