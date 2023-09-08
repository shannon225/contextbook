package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public abstract class AbstractSearchJobData implements SearchJobData {
	private final File diaFile;
	private final PercolatorExecutionData percolatorFiles;
	private final SearchParameters parameters;
	private final String version;
	private String originalDiaFileName=null;

	public AbstractSearchJobData(File diaFile, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version) {
		this.diaFile=diaFile;
		this.percolatorFiles=percolatorFiles;
		this.parameters=parameters;
		this.version=version;
	}
	
	protected File getDiaFile() {
		return diaFile;
	}
	
	@Override
	public String toString() {
		return diaFile.getName();
	}

	@Override
	public StripeFileInterface getDiaFileReader() {
		return StripeFileGenerator.getFile(diaFile, getParameters());
	}
	
	@Override
	public synchronized String getOriginalDiaFileName() {
		if (originalDiaFileName==null) {
			StripeFileInterface file=getDiaFileReader();
			originalDiaFileName=file.getOriginalFileName();
			file.close();
		}
		return originalDiaFileName;
	}
	
	@Override
	public PercolatorExecutionData getPercolatorFiles() {
		// TODO Auto-generated method stub
		return percolatorFiles;
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
		if (!diaFile.exists()) {
			Logger.errorLine("Missing .DIA file: "+diaFile.getName());
			return false;
		}
		if (!percolatorFiles.getInputTSV().exists()) {
			Logger.errorLine("Missing feature file: "+percolatorFiles.getInputTSV().getName());
			return false;
		}
		if (!percolatorFiles.getPeptideOutputFile().exists()) {
			Logger.errorLine("Missing output file: "+percolatorFiles.getPeptideOutputFile().getName());
			return false;
		}
		/*if (!percolatorFiles.getProteinOutputFile().exists()) {
			Logger.errorLine("Missing output file: "+percolatorFiles.getProteinOutputFile().getName());
			return false;
		}*/
		return true;
	}
}