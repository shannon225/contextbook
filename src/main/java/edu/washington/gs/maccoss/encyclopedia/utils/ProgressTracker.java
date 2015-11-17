package edu.washington.gs.maccoss.encyclopedia.utils;

import java.text.DecimalFormat;

public class ProgressTracker {
	DecimalFormat myFormatter = new DecimalFormat("#.##");

	private volatile float progress=0.0f;
	private volatile long lastNotice=0;
	
	public static void main(String[] args) throws InterruptedException {
		ProgressTracker tracker=new ProgressTracker();
		int taskLength=1000;
		float increment=1.0f/taskLength;
		for (int i=0; i<taskLength; i++) {
			tracker.increment(increment);
			Thread.sleep((int)(Math.random()*20));
		}
	}
	
	public void increment(float increment) {
		progress+=increment;
		
		long current=System.currentTimeMillis();
		if (current-lastNotice>500) {
			print();
			lastNotice=current;
		}
	}
	
	public void print() {
		StringBuilder sb=new StringBuilder("\r");
		int length=(int)(progress*100.0f);
		for (int i=0; i<length; i++) {
			sb.append('.');
		}
		sb.append(' ');
		sb.append(myFormatter.format(100.0f*progress));
		sb.append('%');
		System.out.print(sb.toString());
	}
}
