package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import com.sun.istack.Nullable;

import edu.washington.gs.maccoss.encyclopedia.algorithms.OverlapDeconvoluter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class MzmlToDIAConverter implements StripeFileReaderInterface {
	public static final String MZML_EXTENSION=".mzml";

	public static void main(String[] args) throws IOException {
		boolean copy=false;
		
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		//paramMap.put("-filterPeaklists", "true"); 
		
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);
		System.out.println(parameters);
		
		File dir=new File("/Volumes/BriansSSD/bruker/");
		
		File[] files=dir.listFiles(new SimpleFilenameFilter(MZML_EXTENSION));
		for (int i=0; i<files.length; i++) {
			System.out.println((i+1)+" / "+files.length+"\t Copying "+files[i].getName()+"...");

			File f;
			if (copy) {
				f=File.createTempFile(files[i].getName(), MZML_EXTENSION);
				f.deleteOnExit();
				Files.copy(files[i].toPath(), f.toPath(), StandardCopyOption.REPLACE_EXISTING);
			} else {
				f=files[i];
			}
			
			Long time=System.currentTimeMillis();
			File diaFile=new File(files[i].getAbsolutePath().substring(0, files[i].getAbsolutePath().lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
			System.out.println("Converting to "+diaFile.getAbsolutePath());
			
			convertSAX(f, diaFile, parameters, false);

			if (copy) {
				f.delete();
			}
			System.out.println("Total time: "+(System.currentTimeMillis()-time)/1000f+" seconds");
		}
	}
	
	@Override
	public boolean accept(File dir, String name) {
		return name.toLowerCase().endsWith(MzmlToDIAConverter.MZML_EXTENSION);
	}
	
	@Override
	public boolean canTryToReadFile(File f) {
		if (!f.exists()||!f.isFile()||!f.canRead()) return false; 
		return f.getName().toLowerCase().endsWith(MzmlToDIAConverter.MZML_EXTENSION);
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

	static StripeFileInterface convertSAX(File mzMLFile, File diaFile, SearchParameters parameters, boolean isOpenFileInPlace) {
		try {
			Logger.logLine("Indexing "+mzMLFile.getName()+" ...");
			StripeFile stripeFile=new StripeFile(isOpenFileInPlace);
			stripeFile.openFile();

			final BlockingQueue<MSMSBlock> mzmlBlockQueue=new ArrayBlockingQueue<MSMSBlock>(1);
			final MzmlSAXToMSMSProducer producer=new MzmlSAXToMSMSProducer(mzMLFile, 0, mzmlBlockQueue, parameters);

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
				BlockingQueue<MSMSBlock> deconvolutionBlockQueue=new ArrayBlockingQueue<MSMSBlock>(1);
				deconvoluter = new OverlapDeconvoluter(parameters.getFragmentTolerance(), mzmlBlockQueue, deconvolutionBlockQueue);
				retentionTimesByStripe=deconvoluter.getRetentionTimesByStripe();
				ionInjectionTimesByStripe=deconvoluter.getIonInjectionTimesByStripe();
				consumer = new MSMSToDIAConsumer(deconvolutionBlockQueue, stripeFile, parameters);

				Logger.logLine("Converting "+mzMLFile.getName()+" ...");
				deconvoluterThread = new Thread(deconvoluter);
				consumerThread = new Thread(consumer);
			} else {
				deconvoluter = null;
				deconvoluterThread = null;

				consumer = new MSMSToDIAConsumer(mzmlBlockQueue, stripeFile, parameters);

				Logger.logLine("Converting "+mzMLFile.getName()+" ...");
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
				stripeFile.setFileName(mzMLFile.getName(), producer.getMzMLID(), mzMLFile.getAbsolutePath());

				Logger.logLine("Finalizing "+diaFile.getName()+" ...");
				HashMap<Range, WindowData> dutyCycleMap=new HashMap<Range, WindowData>();
				for (Entry<Range, TFloatArrayList> entry : retentionTimesByStripe.entrySet()) {
					Range range=entry.getKey();
					TFloatArrayList rts=entry.getValue();
					float[] deltas=General.firstDerivative(rts.toArray());
					float averageDutyCycle=General.mean(deltas);
					dutyCycleMap.put(range, new WindowData(averageDutyCycle, rts.size()));
				}
				stripeFile.setRanges(dutyCycleMap);

				stripeFile.saveAsFile(diaFile);
				stripeFile.close();
				
				stripeFile=new StripeFile();
				stripeFile.openFile(diaFile);
				Logger.logLine("Finished writing "+diaFile.getName()+"!");

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
	}
}
