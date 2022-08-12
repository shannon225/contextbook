package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import com.sun.istack.Nullable;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors.OverlapDeconvoluter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class DIAProcessor implements StripeFileReaderInterface {
	static final int DEFAULT_QUEUE_CAPACITY = 2;

	
	@Override
	public boolean accept(File dir, String name) {
		return name.toLowerCase().endsWith(StripeFile.DIA_EXTENSION);
	}
	
	@Override
	public boolean canTryToReadFile(File f) {
		if (!f.exists()||!f.isFile()||!f.canRead()) return false; 
		return f.getName().toLowerCase().endsWith(StripeFile.DIA_EXTENSION);
	}
	
	@Override
	public StripeFileInterface readStripeFile(File f, SearchParameters parameters, boolean isOpenFileInPlace) {
		if (canTryToReadFile(f)) {
			String absolutePath=f.getAbsolutePath();
			File diaFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
			return MzmlToDIAConverter.convertSAX(f, diaFile, parameters, isOpenFileInPlace);
		} else {
			throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
		}
	}

	/**
	 * @param inputDIAFile The DIA to convert.
	 * @param outputDIAFile The location where the .DIA file will be saved.
	 * @param parameters Parameters to use during conversion.
	 * @param isOpenFileInPlace TODO: must be true!
	 */
	static StripeFileInterface convertSAX(File inputDIAFile, File outputDIAFile, SearchParameters parameters, boolean isOpenFileInPlace) {
		return convertSAX(inputDIAFile, outputDIAFile, parameters, isOpenFileInPlace, DEFAULT_QUEUE_CAPACITY);
	}

	/**
	 * @param inputDIAFile The input DIA to convert.
	 * @param outputDIAFile The location where the .DIA file will be saved.
	 * @param parameters Parameters to use during conversion.
	 * @param isOpenFileInPlace TODO: must be true!
	 * @param queueCapacity The number of {@link MSMSBlock} kept in the queue(s) between threads. If this many blocks
	 *                      are still pending processing then the processor will block until space is available. Too-low
	 *                      of a setting will lower thread utilization, while too high will use excessive memory.
	 */
	static StripeFileInterface convertSAX(File inputDIAFile, File outputDIAFile, SearchParameters parameters, boolean isOpenFileInPlace, int queueCapacity) {

		/*
		try {
			Logger.logLine("Indexing "+inputDIAFile.getName()+" ...");
			StripeFile stripeFile=new StripeFile(isOpenFileInPlace);
			stripeFile.openFile();

			final BlockingQueue<MSMSBlock> mzmlBlockQueue=new ArrayBlockingQueue<MSMSBlock>(queueCapacity);
			final DIAtoMSMSProducer producer=new DIAtoMSMSProducer(inputDIAFile, 0, mzmlBlockQueue, parameters);

			@Nullable OverlapDeconvoluter deconvoluter;
			MSMSToDIAConsumer consumer;

			final Thread producerThread=new Thread(producer);
			@Nullable final Thread deconvoluterThread;
			Thread consumerThread;

			// will be populated after we join back up. Since we're not looking
			// at it until after the join, we're safe to not have to worry about
			// concurrency.
			HashMap<Range, TFloatArrayList> retentionTimesByStripe=producer.getRetentionTimesByStripe();
			HashMap<Range, TFloatArrayList> ionInjectionTimesByStripe=producer.getIonInjectionTimesByStripe();

			if (parameters.isDeconvoluteOverlappingWindows()) {
				BlockingQueue<MSMSBlock> deconvolutionBlockQueue=new ArrayBlockingQueue<MSMSBlock>(queueCapacity);
				deconvoluter = new OverlapDeconvoluter(parameters.getFragmentTolerance(), mzmlBlockQueue, deconvolutionBlockQueue);
				retentionTimesByStripe=deconvoluter.getRetentionTimesByStripe();
				ionInjectionTimesByStripe=deconvoluter.getIonInjectionTimesByStripe();
				consumer = new MSMSToDIAConsumer(deconvolutionBlockQueue, stripeFile, parameters);

				Logger.logLine("Converting "+inputDIAFile.getName()+" ...");
				deconvoluterThread = new Thread(deconvoluter);
				consumerThread = new Thread(consumer);
			} else {
				deconvoluter = null;
				deconvoluterThread = null;

				consumer = new MSMSToDIAConsumer(mzmlBlockQueue, stripeFile, parameters);

				Logger.logLine("Converting "+inputDIAFile.getName()+" ...");
				consumerThread = new Thread(consumer);

			}

			producerThread.start();
			if (null != deconvoluterThread) {
				deconvoluterThread.start();
			}
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
						if (null != deconvoluterThread) {
							deconvoluterThread.interrupt();
						}
						consumerThread.interrupt();

						if (null != deconvoluterThread) {
							deconvoluterThread.join();
						}
						consumerThread.join();

						throw new EncyclopediaException(producer.getError());
					}

					if (null != deconvoluterThread) {
						deconvoluterThread.join(100L);
						if (deconvoluterThread.isAlive()) {
							alive = true;
						} else if (deconvoluter.hadError()) {
							producerThread.interrupt();
							consumerThread.interrupt();

							producerThread.join();
							consumerThread.join();

							throw new EncyclopediaException(deconvoluter.getError());
						}
					}

					consumerThread.join(100L);
					if (consumerThread.isAlive()) {
						alive = true;
					} else if (consumer.hadError()) {
						producerThread.interrupt();
						if (null != deconvoluterThread) {
							deconvoluterThread.interrupt();
						}

						producerThread.join();
						if (null != deconvoluterThread) {
							deconvoluterThread.join();
						}

						throw new EncyclopediaException(consumer.getError());
					}

					if (!alive) {
						break;
					}
				}

				stripeFile.setSoftwareVersions(producer.getSoftwareAccessionIdToVersion());
				stripeFile.setInstrumentConfiguration(producer.getInstrumentConfigurations());
				if (producer.getStartTime()!=null) {
					// can be missing
					stripeFile.setStartTime(producer.getStartTime());
				}
				stripeFile.setFileName(inputDIAFile.getName(), producer.getMzMLID(), inputDIAFile.getAbsolutePath());

				Logger.logLine("Finalizing "+outputDIAFile.getName()+" ...");
				HashMap<Range, WindowData> dutyCycleMap=new HashMap<Range, WindowData>();
				for (Entry<Range, TFloatArrayList> entry : retentionTimesByStripe.entrySet()) {
					Range range=entry.getKey();
					TFloatArrayList rts=entry.getValue();
					float[] deltas=General.firstDerivative(rts.toArray());
					float averageDutyCycle=General.mean(deltas);
					dutyCycleMap.put(range, new WindowData(averageDutyCycle, rts.size()));
				}
				stripeFile.setRanges(dutyCycleMap);

				stripeFile.saveAsFile(outputDIAFile);
				stripeFile.close();
				
				stripeFile=new StripeFile();
				stripeFile.openFile(outputDIAFile);
				Logger.logLine("Finished writing "+outputDIAFile.getName()+"!");

				return stripeFile;
			} catch (InterruptedException ie) {
				Logger.errorLine("DIA writing interrupted!");
				Logger.errorException(ie);
				return null;
			}
		} catch (IOException ioe) {
			throw new EncyclopediaException("DIA writing IO error!", ioe);
		} catch (SQLException sqle) {
			sqle.printStackTrace();
			throw new EncyclopediaException("DIA writing SQL error!", sqle);
		}
		*/
		return null; //FIXME
	}
}
