package edu.washington.gs.maccoss.encyclopedia.gui.pecan;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.SearchToBLIB;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public class PecanToBLIBJob extends SwingJob {
	private final File blibFile;

	public PecanToBLIBJob(File blibFile, JobProcessor processor) {
		super(processor);
		this.blibFile=blibFile;
	}
	
	@Override
	public String getJobTitle() {
		return "Write "+blibFile.getName();
	}

	@Override
	public void runJob() throws Exception {
		ArrayList<SearchJobData> jobData=new ArrayList<SearchJobData>();
		for (SwingJob job : processor.getQueue()) {
			if (job instanceof PecanJob) {
				jobData.add(((PecanJob)job).getPecanData());
			}
		}

		SearchToBLIB.convert(getProgressIndicator(), jobData, blibFile);
	}
}
