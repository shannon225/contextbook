package edu.washington.gs.maccoss.encyclopedia.jobs;

import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public interface WorkerJob {
	public abstract void runJob(ProgressIndicator progress) throws Exception;
	public abstract String getJobTitle();
}
