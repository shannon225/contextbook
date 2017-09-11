package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public abstract class AbstractSearchJobData implements SearchJobData {

	private final File diaFile;
	private final File featureFile;
	private final File outputFile;
	private final File decoyFile;
	private final SearchParameters parameters;
	private final String version;
	
	public AbstractSearchJobData(File diaFile, File featureFile, File outputFile, File decoyFile, SearchParameters parameters, String version) {
		this.diaFile=diaFile;
		this.featureFile=featureFile;
		this.outputFile=outputFile;
		this.decoyFile=decoyFile;
		this.parameters=parameters;
		this.version=version;
	}

	@Override
	public File getDiaFile() {
		return diaFile;
	}

	@Override
	public File getFeatureFile() {
		return featureFile;
	}

	@Override
	public File getOutputFile() {
		return outputFile;
	}

	@Override
	public File getOutputDecoyFile() {
		return decoyFile;
	}

	@Override
	public SearchParameters getParameters() {
		return parameters;
	}

	@Override
	public String getVersion() {
		return version;
	}

	@Override
	public boolean hasBeenRun() {
		if (!diaFile.exists()) Logger.errorLine("Missing .DIA file: "+diaFile.getName());
		if (!featureFile.exists()) Logger.errorLine("Missing feature file: "+featureFile.getName());
		if (!outputFile.exists()) Logger.errorLine("Missing output file: "+outputFile.getName());
		return diaFile.exists()&&featureFile.exists()&&outputFile.exists();
	}
}