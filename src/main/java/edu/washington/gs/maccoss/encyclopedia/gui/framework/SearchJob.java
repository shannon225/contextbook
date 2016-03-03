package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessor;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SwingJob;

public abstract class SearchJob extends SwingJob {

	protected final SearchJobData data;

	public SearchJob(JobProcessor processor, SearchJobData data) {
		super(processor);
		this.data=data;
	}

	@Override
	public String getJobTitle() {
		return "Read "+data.getDiaFile().getName();
	}

	public SearchJobData getSearchData() {
		return data;
	}

}