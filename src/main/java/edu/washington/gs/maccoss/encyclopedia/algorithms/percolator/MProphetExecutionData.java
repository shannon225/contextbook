package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public class MProphetExecutionData {
	private final File inputTSV;
	private final File fastaFile;
	private final File peptideOutputFile;
	private final File peptideDecoyFile;
	private final SearchParameters parameters;
	public MProphetExecutionData(File inputTSV, File fastaFile, File peptideOutputFile, File peptideDecoyFile, SearchParameters parameters) {
		this.inputTSV = inputTSV;
		this.fastaFile = fastaFile;
		this.peptideOutputFile = peptideOutputFile;
		this.peptideDecoyFile = peptideDecoyFile;
		this.parameters = parameters;
	}
	public File getInputTSV() {
		return inputTSV;
	}
	public File getFastaFile() {
		return fastaFile;
	}
	public File getPeptideOutputFile() {
		return peptideOutputFile;
	}
	public File getPeptideDecoyFile() {
		return peptideDecoyFile;
	}
	public SearchParameters getParameters() {
		return parameters;
	}
}
