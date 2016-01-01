package edu.washington.gs.maccoss.encyclopedia.utils.io;

public class OutputMessage {
	private final String message;
	private final boolean isStdOutput;

	public OutputMessage(String message, boolean isStdOutput) {
		this.message=message;
		this.isStdOutput=isStdOutput;
	}
	
	@Override
	public String toString() {
		return message;
	}
	
	public boolean isStdOutput() {
		return isStdOutput;
	}
	public String getMessage() {
		return message;
	}
}
