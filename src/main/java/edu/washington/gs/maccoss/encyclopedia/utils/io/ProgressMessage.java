package edu.washington.gs.maccoss.encyclopedia.utils.io;

public class ProgressMessage extends OutputMessage {
	private final float progress;

	public ProgressMessage(String message, boolean isStdOutput, float progress) {
		super(message, isStdOutput);
		this.progress=progress;
	}

	public ProgressMessage(String message, float progress) {
		super(message, true);
		this.progress=progress;
	}

	public float getProgress() {
		return progress;
	}
	
	public boolean isError() {
		return progress<0.0f;
	}
	
	public boolean isFinished() {
		return progress>=1.0f;
	}
}
