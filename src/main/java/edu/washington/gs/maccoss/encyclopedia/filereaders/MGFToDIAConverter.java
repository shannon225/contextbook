package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import com.sun.istack.Nullable;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.spectrumprocessors.OverlapDeconvoluter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class MGFToDIAConverter implements StripeFileReaderInterface {
	public static final String MGF_EXTENSION=".mgf";

	public static void main(String[] args) throws IOException {
		boolean copy=false;
		
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);
		
		File dir=new File("/Users/searleb/Downloads/pt_test");
		
		File[] files=dir.listFiles(new SimpleFilenameFilter(MGF_EXTENSION));
		for (int i=0; i<files.length; i++) {
			System.out.println((i+1)+" / "+files.length+"\t Copying "+files[i].getName()+"...");

			File f;
			if (copy) {
				f=File.createTempFile(files[i].getName(), MGF_EXTENSION);
				f.deleteOnExit();
				Files.copy(files[i].toPath(), f.toPath(), StandardCopyOption.REPLACE_EXISTING);
			} else {
				f=files[i];
			}
			
			Long time=System.currentTimeMillis();
			File diaFile=new File(files[i].getAbsolutePath().substring(0, files[i].getAbsolutePath().lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
			System.out.println("Converting to "+diaFile.getAbsolutePath());
			
			convert(f, diaFile, parameters, false);

			if (copy) {
				f.delete();
			}
			System.out.println("Total time: "+(System.currentTimeMillis()-time)/1000f+" seconds");
		}
	}
	
	@Override
	public boolean accept(File dir, String name) {
		return name.toLowerCase().endsWith(MGFToDIAConverter.MGF_EXTENSION);
	}
	
	@Override
	public boolean canTryToReadFile(File f) {
		if (!f.exists()||!f.isFile()||!f.canRead()) return false; 
		return f.getName().toLowerCase().endsWith(MGFToDIAConverter.MGF_EXTENSION);
	}
	
	@Override
	public StripeFileInterface readStripeFile(File f, SearchParameters parameters, boolean isOpenFileInPlace) {
		if (canTryToReadFile(f)) {
			String absolutePath=f.getAbsolutePath();
			File diaFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
			return MGFToDIAConverter.convert(f, diaFile, parameters, isOpenFileInPlace);
		} else {
			throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
		}
	}

	static StripeFileInterface convert(File mgfFile, File diaFile, SearchParameters parameters, boolean isOpenFileInPlace) {
		try {
			Logger.logLine("Indexing "+mgfFile.getName()+" ...");
			StripeFile stripeFile=new StripeFile(isOpenFileInPlace);
			stripeFile.openFile();

			final BlockingQueue<MSMSBlock> mgfBlockQueue=new ArrayBlockingQueue<MSMSBlock>(1);
			final MGFtoMSMSProducer producer=new MGFtoMSMSProducer(mgfFile, mgfBlockQueue);

			Logger.logLine("Converting "+mgfFile.getName()+" ...");
			@Nullable OverlapDeconvoluter deconvoluter;
			MSMSToDIAConsumer consumer=new MSMSToDIAConsumer(mgfBlockQueue, stripeFile, parameters);

			final Thread producerThread=new Thread(producer);
			Thread consumerThread=new Thread(consumer);

			producerThread.start();
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
						consumerThread.interrupt();
						consumerThread.join();

						throw new EncyclopediaException(producer.getError());
					}

					consumerThread.join(100L);
					if (consumerThread.isAlive()) {
						alive = true;
					} else if (consumer.hadError()) {
						producerThread.interrupt();
						producerThread.join();

						throw new EncyclopediaException(consumer.getError());
					}

					if (!alive) {
						break;
					}
				}
				
				stripeFile.setFileName(mgfFile.getName(), mgfFile.getName(), mgfFile.getAbsolutePath());

				Logger.logLine("Finalizing "+diaFile.getName()+" ...");

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
