package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.io.File;
import java.util.ArrayList;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;

public class PecanJobData {
	private Optional<ArrayList<FastaEntry>> targetList;
	private File diaFile;
	private File fastaFile;
	private File featureFile;
	private File outputFile;
	private PecanScoringFactory taskFactory;

	public PecanJobData(Optional<ArrayList<FastaEntry>> targetList, File diaFile, File fastaFile, File featureFile, File outputFile, PecanScoringFactory taskFactory) {
		this.targetList=targetList;
		this.diaFile=diaFile;
		this.fastaFile=fastaFile;
		this.featureFile=featureFile;
		this.outputFile=outputFile;
		this.taskFactory=taskFactory;
	}

	public Optional<ArrayList<FastaEntry>> getTargetList() {
		return targetList;
	}

	public File getDiaFile() {
		return diaFile;
	}

	public File getFastaFile() {
		return fastaFile;
	}

	public File getFeatureFile() {
		return featureFile;
	}

	public File getOutputFile() {
		return outputFile;
	}

	public PecanScoringFactory getTaskFactory() {
		return taskFactory;
	}
}