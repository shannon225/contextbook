package edu.washington.gs.maccoss.encyclopedia.filereaders;

public class WindowData {

	private final float averageDutyCycle;
	private final int numberOfMSMS;

	public WindowData(float averageDutyCycle, int numberOfMSMS) {
		this.averageDutyCycle = averageDutyCycle;
		this.numberOfMSMS = numberOfMSMS;
	}

	public float getAverageDutyCycle() {
		return averageDutyCycle;
	}

	public int getNumberOfMSMS() {
		return numberOfMSMS;
	}
}
