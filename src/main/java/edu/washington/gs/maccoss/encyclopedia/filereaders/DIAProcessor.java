package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors.SpectrumProcessor;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.list.array.TFloatArrayList;

public class DIAProcessor {
	static final int DEFAULT_QUEUE_CAPACITY = 1;
	private final SpectrumProcessor processor;
	private final SearchParameters parameters;
	
	public DIAProcessor(SpectrumProcessor processor, SearchParameters parameters) {
		this.processor = processor;
		this.parameters = parameters;
	}

	public void processStripeFile(ProgressIndicator progress, StripeFileInterface startFile, File resultFile, boolean isOpenFileInPlace) throws IOException, SQLException, DataFormatException {
		StripeFile writeFile=new StripeFile(isOpenFileInPlace);
		writeFile.openFile();

		final BlockingQueue<MSMSBlock> readingBlockQueue=new ArrayBlockingQueue<MSMSBlock>(DEFAULT_QUEUE_CAPACITY);
		final DIAtoMSMSProducer producer=new DIAtoMSMSProducer(startFile, readingBlockQueue);

		final Thread producerThread=new Thread(producer);

		BlockingQueue<MSMSBlock> writingBlockQueue=new ArrayBlockingQueue<MSMSBlock>(DEFAULT_QUEUE_CAPACITY);
		processor.initialize(startFile, parameters, readingBlockQueue, writingBlockQueue);
		
		// not populated yet, but will be at the end
		HashMap<String, String> metadata=producer.getMetadata();
		HashMap<Range, TFloatArrayList> retentionTimesByStripe=processor.getRetentionTimesByStripe();
		HashMap<Range, TFloatArrayList> ionInjectionTimesByStripe=processor.getIonInjectionTimesByStripe();
		
		MSMSToDIAConsumer consumer = new MSMSToDIAConsumer(writingBlockQueue, writeFile, parameters);

		Logger.logLine("Converting "+startFile.getOriginalFileName()+" ...");
		final Thread processorThread = new Thread(processor);
		Thread consumerThread = new Thread(consumer);
		
		producerThread.start();
		processorThread.start();
		consumerThread.start();
		
		try {
			// Spin on threads waiting for execution to finish.
			// If we encounter an error, interrupt the other threads
			// and throw an exception after they die. Otherwise wait
			// until all the threads have finished.
			while (true) {
				boolean alive = false;

				producerThread.join(100L);
				if (producerThread.isAlive()) {
					alive = true;
				} else if (producer.hadError()) {
					processorThread.interrupt();
					consumerThread.interrupt();

					processorThread.join();
					consumerThread.join();

					throw new EncyclopediaException(producer.getError());
				}

				processorThread.join(100L);
				if (processorThread.isAlive()) {
					alive = true;
				} else if (processor.hadError()) {
					producerThread.interrupt();
					consumerThread.interrupt();

					producerThread.join();
					consumerThread.join();

					throw new EncyclopediaException(processor.getError());
				}

				consumerThread.join(100L);
				if (consumerThread.isAlive()) {
					alive = true;
				} else if (consumer.hadError()) {
					producerThread.interrupt();
					processorThread.interrupt();

					producerThread.join();
					processorThread.join();

					throw new EncyclopediaException(consumer.getError());
				}

				if (!alive) {
					break;
				}
			}
			
			writeFile.addMetadata(metadata);

			Logger.logLine("Finalizing "+resultFile.getName()+" ...");
			HashMap<Range, WindowData> dutyCycleMap=new HashMap<Range, WindowData>();
			for (Entry<Range, TFloatArrayList> entry : retentionTimesByStripe.entrySet()) {
				Range range=entry.getKey();
				TFloatArrayList rts=entry.getValue();
				float[] deltas=General.firstDerivative(rts.toArray());
				float averageDutyCycle=General.mean(deltas);
				dutyCycleMap.put(range, new WindowData(averageDutyCycle, rts.size()));
			}
			writeFile.setRanges(dutyCycleMap);

			writeFile.saveAsFile(resultFile);
			writeFile.close();
			Logger.logLine("Finished writing "+resultFile.getName()+"!");

		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		}
	}
}
