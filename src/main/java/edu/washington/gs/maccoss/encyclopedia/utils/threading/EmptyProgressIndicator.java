package edu.washington.gs.maccoss.encyclopedia.utils.threading;

public class EmptyProgressIndicator implements ProgressIndicator {
	volatile private float totalProgress=0.0f;
	
	@Override
	public void update(String message, float totalProgress) {
		System.out.println(((int)(totalProgress*100))+"%\t"+message);
		this.totalProgress=totalProgress;
	}
	
	@Override
	public float getTotalProgress() {
		return totalProgress;
	}
}
