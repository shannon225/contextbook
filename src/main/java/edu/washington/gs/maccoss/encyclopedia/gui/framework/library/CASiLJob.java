package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import edu.washington.gs.maccoss.encyclopedia.CASiL;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;

public class CASiLJob extends SearchJob {
	public CASiLJob(JobProcessor processor, EncyclopediaJobData libraryData) {
		super(processor, libraryData);
	}
	
	@Override
	public void runJob() throws Exception {
		CASiL.runSearch(getProgressIndicator(), getLibraryData());
	}
	
	public EncyclopediaJobData getLibraryData() {
		return (EncyclopediaJobData)getSearchData();
	}
}
