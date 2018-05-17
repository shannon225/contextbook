package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class LineParserConsumer implements Runnable {
	private final BlockingQueue<String> blockingQueue;
	private final LineParserMuscle muscle;

	public LineParserConsumer(BlockingQueue<String> blockingQueue, LineParserMuscle muscle) {
		this.blockingQueue=blockingQueue;
		this.muscle=muscle;
	}
	
	@Override
	public void run() {
		try {
			while (true) {
				String row=blockingQueue.take();
				if (LineParserProducer.POISON_BLOCK==row) break;
				muscle.processRow(row);
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("Table parsing interrupted!");
			Logger.errorException(ie);
		}
	}
}
