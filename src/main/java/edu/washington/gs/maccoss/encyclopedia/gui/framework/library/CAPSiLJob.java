package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import edu.washington.gs.maccoss.encyclopedia.CAPSiL;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;

public class CAPSiLJob extends SearchJob {
	public CAPSiLJob(JobProcessor processor, EncyclopediaJobData libraryData) {
		super(processor, libraryData);
	}
	
	@Override
	public void runJob() throws Exception {
		CAPSiL.runSearch(getProgressIndicator(), getLibraryData());
	}
	
	public EncyclopediaJobData getLibraryData() {
		return (EncyclopediaJobData)getSearchData();
	}
}
