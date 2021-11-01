package edu.washington.gs.maccoss.encyclopedia.gui.framework.scribe;

import edu.washington.gs.maccoss.encyclopedia.Scribe;
import edu.washington.gs.maccoss.encyclopedia.algorithms.scribe.ScribeJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;

public class ScribeJob extends SearchJob {
	public ScribeJob(JobProcessor processor, ScribeJobData pecanData) {
		super(processor, pecanData);
	}
	
	@Override
	public void runJob() throws Exception {
		Scribe.runSearch(getProgressIndicator(), getScribeData());
	}
	
	public ScribeJobData getScribeData() {
		return (ScribeJobData)getSearchData();
	}
}
