package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.util.Map;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
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
		int n = 1; // start at row 1 to account for header line in file;
		           // this way the error message can match line numbering
		try {
			while (true) {
				Map<String, String> row=blockingQueue.take();
				n+=1;
				if (TableParserProducer.POISON_BLOCK==row) break;
				try {
					muscle.processRow(row);
				} catch (Exception e) {
					Logger.errorLine("Exception found parsing table (row="+n+"): "+e.getMessage());
					throw new EncyclopediaException("Error parsing row " + n, e);
				}
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("Table parsing interrupted!");
			Logger.errorException(ie);
		}
	}
}
