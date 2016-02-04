package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

public class SearchJobData {

	private final File diaFile;
	private final File featureFile;
	private final File outputFile;
	private final SearchParameters parameters;
	private final String version;
	
	public SearchJobData(File diaFile, File featureFile, File outputFile, SearchParameters parameters, String version) {
		this.diaFile=diaFile;
		this.featureFile=featureFile;
		this.outputFile=outputFile;
		this.parameters=parameters;
		this.version=version;
	}

	public File getDiaFile() {
		return diaFile;
	}

	public File getFeatureFile() {
		return featureFile;
	}

	public File getOutputFile() {
		return outputFile;
	}

	public SearchParameters getParameters() {
		return parameters;
	}

	public String getVersion() {
		return version;
	}

	public boolean hasBeenRun() {
		return diaFile.exists()&&featureFile.exists()&&outputFile.exists();
	}

}