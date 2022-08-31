package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class SearchToELIBJob extends SwingJob {
	private final File elibFile;
	private final boolean alignBetweenFiles;

	public SearchToELIBJob(File elibFile, boolean alignBetweenFiles, JobProcessor processor) {
		super(processor);
		this.elibFile=elibFile;
		this.alignBetweenFiles=alignBetweenFiles;
	}
	
	@Override
	public String getJobTitle() {
		return "Write Library "+elibFile.getName();
	}

	@Override
	public void runJob() throws Exception {
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (SwingJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}

		Logger.logLine("Found "+jobData.size()+" jobs in the queue to combine...");
		SearchToBLIB.convert(getProgressIndicator(), jobData, elibFile, false, alignBetweenFiles);
	}
}
