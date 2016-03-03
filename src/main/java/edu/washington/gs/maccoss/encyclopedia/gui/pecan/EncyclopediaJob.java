package edu.washington.gs.maccoss.encyclopedia.gui.pecan;

import edu.washington.gs.maccoss.encyclopedia.Encyclopedia;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;

public class EncyclopediaJob extends SearchJob {
	public EncyclopediaJob(JobProcessor processor, EncyclopediaJobData libraryData) {
		super(processor, libraryData);
	}
	
	@Override
	public void runJob() throws Exception {
		Encyclopedia.runSearch(getProgressIndicator(), getLibraryData());
	}
	
	public EncyclopediaJobData getLibraryData() {
		return (EncyclopediaJobData)getSearchData();
	}
}
