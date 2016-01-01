package edu.washington.gs.maccoss.encyclopedia.utils;

import java.util.ArrayList;

public class Logger {
	private static final ArrayList<LogRecorder> recorders=new ArrayList<LogRecorder>();
	public static void addRecorder(LogRecorder recorder) {
		recorders.add(recorder);
	}
	
	public static void log(String s) {
		System.out.print(s);
		for (LogRecorder recorder : recorders) {
			recorder.log(s);
		}
	}
	public static void logLine(String s) {
		System.out.println(s);
		for (LogRecorder recorder : recorders) {
			recorder.logLine(s);
		}
	}

	public static void errorLine(String s) {
		System.err.println(s);
		for (LogRecorder recorder : recorders) {
			recorder.errorLine(s);
		}
	}

	public static void logException(Throwable e) {
		System.out.println(e);
		for (StackTraceElement ste : e.getStackTrace()) {
			System.out.println("\t"+ste.toString());
		}
		for (LogRecorder recorder : recorders) {
			recorder.logException(e);
		}
	}

	public static void errorException(Throwable e) {
		System.err.println(e);
		for (StackTraceElement ste : e.getStackTrace()) {
			System.err.println("\t"+ste.toString());
		}
		for (LogRecorder recorder : recorders) {
			recorder.errorException(e);
		}
	}
}
