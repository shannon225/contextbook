package edu.washington.gs.maccoss.encyclopedia.gui.framework.library;

import edu.washington.gs.maccoss.encyclopedia.Thesaurus;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;

public class ThesaurusJob extends SearchJob {
	public ThesaurusJob(JobProcessor processor, ThesaurusJobData libraryData) {
		super(processor, libraryData);
	}
	
	@Override
	public void runJob() throws Exception {
		Thesaurus.runSearch(getProgressIndicator(), getLibraryData());
	}
	
	public ThesaurusJobData getLibraryData() {
		return (ThesaurusJobData)getSearchData();
	}
}
