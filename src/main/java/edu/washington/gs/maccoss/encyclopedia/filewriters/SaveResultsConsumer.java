package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class SaveResultsConsumer implements PeptideScoringResultsConsumer {
	// only written to by one thread
	private final ArrayList<PeptideScoringResult> savedResults=new ArrayList<PeptideScoringResult>();
	private final BlockingQueue<PeptideScoringResult> resultsQueue;

	public SaveResultsConsumer(BlockingQueue<PeptideScoringResult> resultsQueue) {
		this.resultsQueue=resultsQueue;
	}
	
	public ArrayList<PeptideScoringResult> getSavedResults() {
		return new ArrayList<PeptideScoringResult>(savedResults);
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
	public BlockingQueue<PeptideScoringResult> getResultsQueue() {
		return resultsQueue;
	}

	@Override
	public void run() {
		try {
			while (true) {
				PeptideScoringResult result=resultsQueue.take();
				if (PeptideScoringResult.POISON_RESULT==result) break;
				savedResults.add(result);
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA processing interrupted!");
			Logger.errorException(ie);
		}

	}
}
