package edu.washington.gs.maccoss.encyclopedia.gui.general;

import java.util.List;

import javax.swing.SwingWorker;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.io.ProgressMessage;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public abstract class SwingJob extends SwingWorker<Nothing, ProgressMessage> {

	protected volatile String message="";
	protected volatile float progress=0.0f;
	protected final JobProcessor processor;

	public SwingJob(JobProcessor processor) {
		super();
		this.processor=processor;
	}

	public abstract void runJob() throws Exception;

	public abstract String getJobTitle();

	@Override
	protected Nothing doInBackground() {
		try {
			runJob();
		} catch (Exception e) {
			publish(new ProgressMessage("Encountered Fatal Error!", -1.0f));
			progress=-1.0f;
			message="Encountered Fatal Error!";
			Logger.errorException(e);
		}
		return Nothing.NOTHING;
	}

	@Override
	protected void done() {
		if (progress>0) { // not in error
			progress=1.0f;
		}
	}

	@Override
	protected void process(List<ProgressMessage> chunks) {
		for (ProgressMessage p : chunks) {
			if (progress<=p.getProgress()) {
				progress=p.getProgress();
				message=p.getMessage();
			}
		}
	}

	public ProgressIndicator getProgressIndicator() {
		final ProgressIndicator indicator=new ProgressIndicator() {
			volatile private float totalProgress=0.0f;
			@Override
			public void update(String message, float totalProgress) {
				publish(new ProgressMessage(message, totalProgress));
				this.totalProgress=totalProgress;
			}
			@Override
			public void update(String message) {
				publish(new ProgressMessage(message, totalProgress));
			}
			@Override
			public float getTotalProgress() {
				return totalProgress;
			}
		};
		return indicator;
	}

	public ProgressMessage getProgressMessage() {
		synchronized (this) {
			return new ProgressMessage(message, progress);
		}
	}

	public float getTotalProgress() {
		return progress;
	}

	public String getMessage() {
		return message;
	}

	public boolean isFinished() {
		return getProgressMessage().isFinished();
	}

	public boolean isError() {
		return getProgressMessage().isError();
	}

}