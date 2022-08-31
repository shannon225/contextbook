package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.util.List;

import edu.washington.gs.maccoss.encyclopedia.jobs.WorkerJob;
import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressMessage;

public class SwingWorkerJob extends SwingJob {
	private final WorkerJob job;
	

	public SwingWorkerJob(JobProcessor processor, WorkerJob job) {
		super(processor);
		this.job=job;
	}

	@Override
	public void runJob() throws Exception {
		job.runJob(getProgressIndicator());
	}

	@Override
	public String getJobTitle() {
		return job.getJobTitle();
	}

	public WorkerJob getJob() {
		return job;
	}

	@Override
	protected void done() {
		super.done();
		processor.fireJobUpdated(getJob());
	}

	@Override
	protected void process(List<ProgressMessage> chunks) {
		boolean needsToUpdate=false;
		for (ProgressMessage p : chunks) {
			if (progress<=p.getProgress()) {
				needsToUpdate=true;
			}
		}
		
		super.process(chunks);
		
		if (needsToUpdate) {
			processor.fireJobUpdated(getJob());
		}
	}
}
