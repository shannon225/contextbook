package edu.washington.gs.maccoss.encyclopedia.jobs;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class SearchToBLIBJob implements WorkerJob {
	private final File blibFile;
	private final boolean alignBetweenFiles;
	JobProcessor processor;

	public SearchToBLIBJob(File blibFile, boolean alignBetweenFiles, JobProcessor processor) {
		this.processor=processor;
		this.blibFile=blibFile;
		this.alignBetweenFiles=alignBetweenFiles;
	}
	
	@Override
	public String getJobTitle() {
		return "Write BLIB "+blibFile.getName();
	}

	@Override
	public void runJob(ProgressIndicator progress) throws Exception {
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (WorkerJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}

		Logger.logLine("Found "+jobData.size()+" jobs in the queue to combine...");
		SearchToBLIB.convert(progress, jobData, blibFile, true, alignBetweenFiles);
	}
}
