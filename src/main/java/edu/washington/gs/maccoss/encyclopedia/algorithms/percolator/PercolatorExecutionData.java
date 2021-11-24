package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public class PercolatorExecutionData {
	private final File inputTSV;
	private final File fastaFile;
	private final File peptideOutputFile;
	private final File peptideDecoyFile;
	private final File proteinOutputFile;
	private final File proteinDecoyFile;
	private final SearchParameters parameters;
	private final boolean useMinMax;
	private String percolatorExecutableVersion;

	public PercolatorExecutionData(File inputTSV, File fastaFile, File peptideOutputFile, File peptideDecoyFile, File proteinOutputFile, File proteinDecoyFile, SearchParameters parameters) {
		this.inputTSV=inputTSV;
		this.fastaFile=fastaFile;
		this.peptideOutputFile=peptideOutputFile;
		this.peptideDecoyFile=peptideDecoyFile;
		this.proteinOutputFile=proteinOutputFile;
		this.proteinDecoyFile=proteinDecoyFile;
		this.parameters=parameters;
		this.useMinMax=true;
	}
	
	
	
	public PercolatorExecutionData(File inputTSV, File fastaFile, File peptideOutputFile, File peptideDecoyFile, File proteinOutputFile, File proteinDecoyFile, SearchParameters parameters, boolean useMinMax) {
		this.inputTSV=inputTSV;
		this.fastaFile=fastaFile;
		this.peptideOutputFile=peptideOutputFile;
		this.peptideDecoyFile=peptideDecoyFile;
		this.proteinOutputFile=proteinOutputFile;
		this.proteinDecoyFile=proteinDecoyFile;
		this.parameters=parameters;
		this.useMinMax=useMinMax;
	}

	public PercolatorExecutionData getDDAVersion() {
		return new PercolatorExecutionData(inputTSV, fastaFile, peptideOutputFile, peptideDecoyFile, proteinOutputFile, proteinDecoyFile, parameters, false);
	}
	
	public boolean hasDataAvailable() {
		if (!peptideOutputFile.exists()||!peptideOutputFile.canRead()) return false;
		if (!peptideDecoyFile.exists()||!peptideDecoyFile.canRead()) return false;
		return true;
	}
	
	public File getFastaFile() {
		return fastaFile;
	}

	public File getInputTSV() {
		return inputTSV;
	}

	public File getWeightsFile(int round) {
		return new File(getPeptideOutputFile().getAbsolutePath()+"."+round+".weights");
	}

	public File getModelFile() {
		return new File(getPeptideOutputFile().getAbsolutePath()+".model");
	}

	/**
	 * This method is {@code protected} to allow access for testing.
	 *
	 * @param percolatorExecutableVersion Canonical version of Percolator parsed directly when running the actual executable
	 */
	protected void setPercolatorExecutableVersion(String percolatorExecutableVersion) {
		this.percolatorExecutableVersion = percolatorExecutableVersion;
	}

	/**
	 * @return Canonical version of Percolator parsed directly when running the actual executable
	 */
	public Optional<String> getPercolatorExecutableVersion() {
		return Optional.ofNullable(this.percolatorExecutableVersion);
	}

	public File getPeptideOutputFile() {
		return peptideOutputFile;
	}

	public File getPeptideDecoyFile() {
		return peptideDecoyFile;
	}

	public File getProteinOutputFile() {
		return proteinOutputFile;
	}

	public File getProteinDecoyFile() {
		return proteinDecoyFile;
	}
	public SearchParameters getParameters() {
		return parameters;
	}
	
	public boolean isUseMinMax() {
		return useMinMax;
	}
	
	public Optional<File> getPercolatorModelFile() {
		return parameters.getPercolatorModelFile();
	}
}
