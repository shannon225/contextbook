package edu.washington.gs.maccoss.encyclopedia.utils.io;

import java.io.File;
import java.io.InputStream;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class LineParser {
	public static void parseFile(InputStream is, LineParserMuscle muscle) {

		BlockingQueue<String> blockingQueue=new LinkedBlockingQueue<String>();
		LineParserProducer producer=new LineParserProducer(blockingQueue, is, 1);
		LineParserConsumer consumer=new LineParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
			
			muscle.cleanup();
			
		} catch (InterruptedException ie) {
			Logger.errorLine("LineParser.parseFile("+muscle.getClass().getName()+") reading interrupted!");
			Logger.errorException(ie);
		}
	}

	public static void parseFile(File f, LineParserMuscle muscle) {

		BlockingQueue<String> blockingQueue=new LinkedBlockingQueue<String>();
		LineParserProducer producer=new LineParserProducer(blockingQueue, f, 1);
		LineParserConsumer consumer=new LineParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
			
			muscle.cleanup();
			
		} catch (InterruptedException ie) {
			Logger.errorLine("LineParser.parseFile("+muscle.getClass().getName()+") reading interrupted!");
			Logger.errorException(ie);
		}
	}
}
