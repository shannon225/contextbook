package edu.washington.gs.maccoss.encyclopedia.gui.framework;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.jobs.WorkerJob;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;

public abstract class SearchJob implements WorkerJob, XMLObject {

	protected final SearchJobData data;

	public SearchJob(SearchJobData data) {
		this.data=data;
	}

	@Override
	public String getJobTitle() {
		return "Read "+data;//.getDiaFileReader().getOriginalFileName();
	}

	public SearchJobData getSearchData() {
		return data;
	}

}