package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class SaveResultsConsumer implements PeptideScoringResultsConsumer {
	// only written to by one thread
	private final ArrayList<AbstractScoringResult> savedResults=new ArrayList<AbstractScoringResult>();
	private final BlockingQueue<AbstractScoringResult> resultsQueue;

	public SaveResultsConsumer(BlockingQueue<AbstractScoringResult> resultsQueue) {
		this.resultsQueue=resultsQueue;
	}
	
	public ArrayList<AbstractScoringResult> getSavedResults() {
		return new ArrayList<AbstractScoringResult>(savedResults);
	}

	@Override
	public void close() {
	}

	@Override
	public int getNumberProcessed() {
		// may be out of date, but at least won't throw concurrency errors
		return savedResults.size();
	}

	@Override
	public BlockingQueue<AbstractScoringResult> getResultsQueue() {
		return resultsQueue;
	}

	@Override
	public void run() {
		try {
			while (true) {
				AbstractScoringResult result=resultsQueue.take();
				if (AbstractScoringResult.POISON_RESULT==result) break;
				savedResults.add(result);
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA processing interrupted!");
			Logger.errorException(ie);
		}

	}
}
