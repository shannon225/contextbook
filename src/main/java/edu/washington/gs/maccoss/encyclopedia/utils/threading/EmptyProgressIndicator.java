package edu.washington.gs.maccoss.encyclopedia.utils.threading;

public class EmptyProgressIndicator implements ProgressIndicator {
	volatile private float totalProgress=0.0f;
	private final boolean print;

	public EmptyProgressIndicator() {
		this.print=false;
	}
	
	public EmptyProgressIndicator(boolean print) {
		this.print=print;
	}

	@Override
	public void update(String message, float totalProgress) {
		if (print) {
			System.out.println(((int)(totalProgress*100))+"%\t"+message);
		}
		this.totalProgress=totalProgress;
	}
	
	@Override
	public void update(String message) {
		if (print) {
			System.out.println(((int)(totalProgress*100))+"%\t"+message);
		}
	}
	
	@Override
	public float getTotalProgress() {
		return totalProgress;
	}
}
