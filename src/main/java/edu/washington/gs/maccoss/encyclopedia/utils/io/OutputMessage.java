package edu.washington.gs.maccoss.encyclopedia.utils.io;

public class OutputMessage {
	public final String message;
	public final boolean isStdOutput;

	public OutputMessage(String message, boolean isStdOutput) {
		this.message=message;
		this.isStdOutput=isStdOutput;
	}
	
	@Override
	public String toString() {
		return message;
	}
}
