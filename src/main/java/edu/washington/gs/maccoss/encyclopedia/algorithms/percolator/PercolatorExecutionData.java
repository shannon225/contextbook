package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.io.File;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.io.Version;

public class PercolatorExecutionData {
	private final File inputTSV;
	private final File fastaFile;
	private final File peptideOutputFile;
	private final File peptideDecoyFile;
	private final File proteinOutputFile;
	private final File proteinDecoyFile;
	private final SearchParameters parameters;
	private Version percolatorExecutableVersion;

	public PercolatorExecutionData(File inputTSV, File fastaFile, File peptideOutputFile, File peptideDecoyFile, File proteinOutputFile, File proteinDecoyFile, SearchParameters parameters) {
		this.inputTSV=inputTSV;
		this.fastaFile=fastaFile;
		this.peptideOutputFile=peptideOutputFile;
		this.peptideDecoyFile=peptideDecoyFile;
		this.proteinOutputFile=proteinOutputFile;
		this.proteinDecoyFile=proteinDecoyFile;
		this.parameters=parameters;
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

	/**
	 * @param percolatorExecutableVersion Canonical version of Percolator parsed directly when running the actual executable
	 */
	void setPercolatorExecutableVersion(Version percolatorExecutableVersion) {
		this.percolatorExecutableVersion = percolatorExecutableVersion;
	}

	/**
	 * @return Canonical version of Percolator parsed directly when running the actual executable
	 */
	public Optional<Version> getPercolatorExecutableVersion() {
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
}
