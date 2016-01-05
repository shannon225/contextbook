package edu.washington.gs.maccoss.encyclopedia.gui.pecan;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.PecanToBLIB;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

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
		ArrayList<PecanJobData> jobData=new ArrayList<PecanJobData>();
		for (SwingJob job : processor.getQueue()) {
			System.out.println(job.getTotalProgress()+", "+job.getJobTitle()+"\t"+job.isFinished());
			if (job.isFinished()&&job instanceof PecanJob) {
				jobData.add(((PecanJob)job).getPecanData());
			}
		}

		PecanToBLIB.convert(new EmptyProgressIndicator(), jobData, blibFile);
	}
}
