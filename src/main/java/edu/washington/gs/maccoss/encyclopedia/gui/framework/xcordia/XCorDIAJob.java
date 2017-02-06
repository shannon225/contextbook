package edu.washington.gs.maccoss.encyclopedia.gui.framework.xcordia;

import edu.washington.gs.maccoss.encyclopedia.XCorDIA;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;

public class XCorDIAJob extends SearchJob {
	public XCorDIAJob(JobProcessor processor, XCorDIAJobData pecanData) {
		super(processor, pecanData);
	}
	
	@Override
	public void runJob() throws Exception {
		XCorDIA.runPie(getProgressIndicator(), getXCorDIAData());
	}
	
	public XCorDIAJobData getXCorDIAData() {
		return (XCorDIAJobData)getSearchData();
	}
}
