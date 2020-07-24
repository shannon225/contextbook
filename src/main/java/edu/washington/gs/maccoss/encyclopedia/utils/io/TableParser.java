package edu.washington.gs.maccoss.encyclopedia.utils.io;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

import java.io.File;
import java.util.Map;
import java.util.concurrent.*;

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

		parseTable(muscle, producer);
	}

	/**
	 * To support unit testing. Will ALWAYS use a single consumer.
	 */
	static void parseTable(TableParserMuscle muscle, TableParserProducer producer) {
		final TableParserConsumer consumer = new TableParserConsumer(producer.getQueue(), muscle);

		// Set up a thread pool for the multithreaded processing
		final ForkJoinPool executor = new ForkJoinPool(2);
		try {
			// Use futures to execute code on the pool
			// and report back any exceptions.

			final CompletableFuture<Void> producerFuture = CompletableFuture.runAsync(producer, executor);
			final CompletableFuture<Void> consumerFuture = CompletableFuture.runAsync(consumer, executor);

			try {
				try {
					producerFuture.get();
				} catch (ExecutionException e) {
					throw new EncyclopediaException("Error reading tabular file", e);
				}
				try {
					consumerFuture.get();
				} catch (ExecutionException e) {
					throw new EncyclopediaException("Error parsing tabular file (" + muscle.getClass().getName() + ")", e);
				}
			} catch (InterruptedException e) {
				throw new EncyclopediaException("Interrupted while parsing tabular file (" + muscle.getClass().getName() + ")", e);
			}
		} finally {
			executor.shutdownNow();

			// Clean up once the thread pool has been shut down
			muscle.cleanup();
		}
	}
}
