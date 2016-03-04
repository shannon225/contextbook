package edu.washington.gs.maccoss.encyclopedia.gui.framework.pecan;

import edu.washington.gs.maccoss.encyclopedia.Pecanpie;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;

public class PecanJob extends SearchJob {
	public PecanJob(JobProcessor processor, PecanJobData pecanData) {
		super(processor, pecanData);
	}
	
	@Override
	public void runJob() throws Exception {
		Pecanpie.runPie(getProgressIndicator(), getPecanData());
	}
	
	public PecanJobData getPecanData() {
		return (PecanJobData)getSearchData();
	}
}
