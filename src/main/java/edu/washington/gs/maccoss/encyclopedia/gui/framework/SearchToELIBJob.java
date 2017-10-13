package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public class SearchToELIBJob extends SwingJob {
	private final File elibFile;
	private final File fastaFile;
	private final boolean alignBetweenFiles;

	public SearchToELIBJob(File elibFile, File fastaFile, boolean alignBetweenFiles, JobProcessor processor) {
		super(processor);
		this.elibFile=elibFile;
		this.fastaFile=fastaFile;
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

		SearchToBLIB.convert(getProgressIndicator(), jobData, fastaFile, elibFile, false, alignBetweenFiles);
	}
}
