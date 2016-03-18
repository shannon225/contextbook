package edu.washington.gs.maccoss.encyclopedia.utils.threading;

public class SubProgressIndicator implements ProgressIndicator {
	private final ProgressIndicator progress;
	private final float start;
	private final float percentage;
	public SubProgressIndicator(ProgressIndicator progress, float percentage) {
		this.progress=progress;
		this.start=progress.getTotalProgress();
		this.percentage=percentage;
	}

	@Override
	public void update(String message, float totalProgress) {
		if (totalProgress>1.0f) totalProgress=1.0f;
		progress.update(message, start+totalProgress*percentage);
	}
	
	@Override
	public void update(String message) {
		progress.update(message);
	}	
	
	@Override
	public float getTotalProgress() {
		return progress.getTotalProgress();
	}
	
	public float getPercentage() {
		return percentage;
	}
}
