package edu.washington.gs.maccoss.encyclopedia.utils.threading;

public class EmptyProgressIndicator implements ProgressIndicator {
	@Override
	public void update(String message, float totalProgress) {
		System.out.println(((int)(totalProgress*100))+"%\t"+message);
	}
}
