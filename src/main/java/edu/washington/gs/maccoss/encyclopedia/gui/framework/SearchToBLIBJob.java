package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public class SearchToBLIBJob extends SwingJob {
	private final File blibFile;
	private final boolean alignBetweenFiles;

	public SearchToBLIBJob(File blibFile, boolean alignBetweenFiles, JobProcessor processor) {
		super(processor);
		this.blibFile=blibFile;
		this.alignBetweenFiles=alignBetweenFiles;
	}
	
	@Override
	public String getJobTitle() {
		return "Write BLIB "+blibFile.getName();
	}

	@Override
	public void runJob() throws Exception {
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (SwingJob job : processor.getQueue()) {
			if (job instanceof SearchJob) {
				jobData.add(((SearchJob)job).getSearchData());
			}
		}

		SearchToBLIB.convert(getProgressIndicator(), jobData, blibFile, true, alignBetweenFiles);
	}
}
