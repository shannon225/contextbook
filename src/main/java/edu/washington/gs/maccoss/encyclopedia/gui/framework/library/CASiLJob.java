package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import edu.washington.gs.maccoss.encyclopedia.Thesaurus;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.CASiLJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;

public class CASiLJob extends SearchJob {
	public CASiLJob(JobProcessor processor, CASiLJobData libraryData) {
		super(processor, libraryData);
	}
	
	@Override
	public void runJob() throws Exception {
		Thesaurus.runSearch(getProgressIndicator(), getLibraryData());
	}
	
	public CASiLJobData getLibraryData() {
		return (CASiLJobData)getSearchData();
	}
}
