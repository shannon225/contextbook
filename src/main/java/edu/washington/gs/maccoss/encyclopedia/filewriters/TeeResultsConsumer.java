package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class TeeResultsConsumer implements PeptideScoringResultsConsumer {

	private final BlockingQueue<PeptideScoringResult> resultsQueue;
	private final PeptideScoringResultsConsumer[] consumers;
	
	public TeeResultsConsumer(BlockingQueue<PeptideScoringResult> resultsQueue, PeptideScoringResultsConsumer... consumers) {
		this.resultsQueue=resultsQueue;
		this.consumers=consumers;
	}
	
	@Override
	public BlockingQueue<PeptideScoringResult> getResultsQueue() {
		return resultsQueue;
	}
	
	@Override
	public void close() {
		for (PeptideScoringResultsConsumer consumer : consumers) {
			consumer.close();
		}
	}
	
	@Override
	public int getNumberProcessed() {
		int numProcessed=Integer.MAX_VALUE;
		for (PeptideScoringResultsConsumer consumer : consumers) {
			int processed=consumer.getNumberProcessed();
			if (processed<numProcessed) {
				numProcessed=processed;
			}
		}
		return numProcessed;
	}
	
	@Override
	public void run() {
		try {
			while (true) {
				PeptideScoringResult result=resultsQueue.take();
				for (PeptideScoringResultsConsumer consumer : consumers) {
					consumer.getResultsQueue().add(result);
				}
				
				// make sure you give the poison result to the children first
				if (PeptideScoringResult.POISON_RESULT==result) break;
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA processing interrupted!");
			Logger.errorException(ie);
		}
		
	}
}
