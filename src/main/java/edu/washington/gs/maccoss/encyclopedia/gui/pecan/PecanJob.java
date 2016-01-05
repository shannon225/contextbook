package edu.washington.gs.maccoss.encyclopedia.gui.pecan;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public class PecanJob extends SwingJob {
	private final PecanJobData pecanData;

	public PecanJob(JobProcessor processor, PecanJobData pecanData) {
		super(processor);
		this.pecanData=pecanData;
	}
	@Override
	public void runJob() throws Exception {
		Pecanpie.runPie(getProgressIndicator(), pecanData);
	}
	@Override
	public String getJobTitle() {
		return "Read "+pecanData.getDiaFile().getName();
	}
	
	public PecanJobData getPecanData() {
		return pecanData;
	}
}
