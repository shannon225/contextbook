package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.util.Map;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class TableParserConsumer implements Runnable {
	private final BlockingQueue<Map<String, String>> blockingQueue;
	private final TableParserMuscle muscle;

	public TableParserConsumer(BlockingQueue<Map<String, String>> blockingQueue, TableParserMuscle muscle) {
		this.blockingQueue=blockingQueue;
		this.muscle=muscle;
	}
	
	@Override
	public void run() {
		try {
			while (true) {
				Map<String, String> row=blockingQueue.take();
				if (TableParserProducer.POISON_BLOCK==row) break;
				muscle.processRow(row);
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("Table parsing interrupted!");
			Logger.errorException(ie);
		}
	}
}
