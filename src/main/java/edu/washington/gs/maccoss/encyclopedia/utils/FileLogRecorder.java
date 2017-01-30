package edu.washington.gs.maccoss.encyclopedia.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class FileLogRecorder implements LogRecorder {
	private final SimpleDateFormat format=new SimpleDateFormat("[HH:mm:ss] ");
	private final PrintWriter writer;
	
	private long lastFlush=0;
	private long flushInterval=1000;
	
	public FileLogRecorder(File f) throws FileNotFoundException, UnsupportedEncodingException {
		writer=new PrintWriter(f, "UTF-8");
	}

	private void maybeFlush() {
		if (System.currentTimeMillis()-lastFlush>flushInterval) {
			writer.flush();
			lastFlush=System.currentTimeMillis();
		}
	}
	
	public void log(String s) {
		writer.print(s);
		maybeFlush();
	}
	public void logLine(String s) {
		writer.println(format.format(new Date())+s);
		maybeFlush();
	}
	public void timelessLogLine(String s) {
		writer.println(s);
		maybeFlush();
	}

	public void errorLine(String s) {
		writer.println("ERROR "+format.format(new Date())+s);
		maybeFlush();
	}

	public void logException(Throwable e) {
		writer.println(format.format(new Date())+e);
		for (StackTraceElement ste : e.getStackTrace()) {
			writer.println("\t"+ste.toString());
		}
		maybeFlush();
	}

	public void errorException(Throwable e) {
		writer.println("ERROR "+format.format(new Date())+e);
		for (StackTraceElement ste : e.getStackTrace()) {
			writer.println("ERROR \t"+ste.toString());
		}
		maybeFlush();
	}
	
	public void close() {
		writer.flush();
		writer.close();
	}
}
