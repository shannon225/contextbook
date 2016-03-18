package edu.washington.gs.maccoss.encyclopedia.utils.threading;

public interface ProgressIndicator {
	public void update(String message);
	public void update(String message, float totalProgress);
	public float getTotalProgress();
}
