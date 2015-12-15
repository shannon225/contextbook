package edu.washington.gs.maccoss.encyclopedia.utils;

public class Logger {
	public static void log(String s) {
		System.out.print(s);
	}
	public static void logLine(String s) {
		System.out.println(s);
	}

	public static void errorLine(String s) {
		System.err.println(s);
	}

	public static void logException(Throwable e) {
		System.out.println(e);
		for (StackTraceElement ste : e.getStackTrace()) {
			System.out.println("\t"+ste.toString());
		}
	}

	public static void errorException(Throwable e) {
		System.err.println(e);
		for (StackTraceElement ste : e.getStackTrace()) {
			System.err.println("\t"+ste.toString());
		}
	}
}
