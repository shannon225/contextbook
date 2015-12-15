package edu.washington.gs.maccoss.encyclopedia.filewriters;

public interface PeptideScoringResultsConsumer extends Runnable {
	public int getNumberProcessed();
	public void close();
}
