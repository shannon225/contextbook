package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;

public interface PeptideScoringResultsConsumer extends Runnable {
	public int getNumberProcessed();
	public void close();
	public BlockingQueue<AbstractScoringResult> getResultsQueue();
}
