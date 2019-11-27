package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.File;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class TableParser {
	public static void parseCSV(File f, TableParserMuscle muscle) {
		parseTable(f, muscle, ",");
	}

	public static void parseTSV(File f, TableParserMuscle muscle) {
		parseTable(f, muscle, "\t");
	}

	public static void parseSSV(File f, TableParserMuscle muscle) {
		parseTable(f, muscle, " ");
	}

	public static void parseTable(File f, TableParserMuscle muscle, String token) {

		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, f, token, 1);
		TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
			
			muscle.cleanup();
			
		} catch (InterruptedException ie) {
			Logger.errorLine("TableParser.parseTable("+muscle.getClass().getName()+") reading interrupted!");
			Logger.errorException(ie);
		}
	}
}
