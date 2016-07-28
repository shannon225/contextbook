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

import edu.washington.gs.maccoss.encyclopedia.algorithms.OverlapDeconvoluter;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

public class MzmlToDIAConverter {
	public static final String MZML_EXTENSION=".mzml";

	public static void main(String[] args) throws IOException {
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);
		System.out.println(parameters);
		
		File dir=new File("/Volumes/WorkingDisk/yeast_curve/");
		
		File[] files=dir.listFiles(new SimpleFilenameFilter(MZML_EXTENSION));
		for (int i=0; i<files.length; i++) {
			System.out.println((i+1)+" / "+files.length+"\t Copying "+files[i].getName()+"...");

			File f=File.createTempFile(files[i].getName(), MZML_EXTENSION);
			f.deleteOnExit();
			Files.copy(files[i].toPath(), f.toPath(), StandardCopyOption.REPLACE_EXISTING);
			
			Long time=System.currentTimeMillis();
			File diaFile=new File(files[i].getAbsolutePath().substring(0, files[i].getAbsolutePath().lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
			System.out.println("Converting to "+diaFile.getAbsolutePath());
			
			convert(f, diaFile, parameters);
			
			f.delete();
			System.out.println("Total time: "+(System.currentTimeMillis()-time)/1000f+" seconds");
		}
	}

	public static StripeFileInterface getFile(File f, SearchParameters parameters) {
		if (!f.exists()||!f.canRead()) {
			throw new EncyclopediaException("Can't read file "+f.getAbsolutePath());
		}
		
		// first try to read if .DIA
		if (f.getName().toLowerCase().endsWith(StripeFile.DIA_EXTENSION)) {
			return openDIAFile(f);
		}
		
		// then try to change name to .DIA and read
		String absolutePath=f.getAbsolutePath();
		File diaFile=new File(absolutePath.substring(0, absolutePath.lastIndexOf('.'))+StripeFile.DIA_EXTENSION);
		if (diaFile.exists()&&diaFile.canRead()) {
			return openDIAFile(diaFile);
		}
		
		// otherwise check for MZML and convert
		if (f.getName().toLowerCase().endsWith(MZML_EXTENSION)) {
			return convert(f, diaFile, parameters);
		} else {
			throw new EncyclopediaException("Can't read file type "+f.getAbsolutePath());
		}
	}

	public static StripeFileInterface openDIAFile(File f) {
		try {
			StripeFileInterface stripefile=new StripeFile();
			stripefile.openFile(f);
			return stripefile;
		} catch (IOException ioe) {
			Logger.errorLine("Unexpected exception reading DIA file: "+f.getName());
			Logger.errorException(ioe);
			throw new EncyclopediaException("Error reading DIA file!", ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Unexpected exception reading DIA file: "+f.getName());
			Logger.logException(sqle);
			throw new EncyclopediaException("Error reading DIA file!", sqle);
		}
	}

	static StripeFileInterface convert(File mzMLFile, File diaFile, SearchParameters parameters) {
		try {
			Logger.logLine("Indexing "+mzMLFile.getName()+" ...");
			StripeFile stripeFile=new StripeFile();
			stripeFile.openFile();

			MzMLUnmarshaller unmarshaller=new MzMLUnmarshaller(mzMLFile);
			stripeFile.setFileName(mzMLFile.getName(), unmarshaller.getMzMLId(), mzMLFile.getAbsolutePath());

			BlockingQueue<MzmlBlock> mzmlBlockQueue=new ArrayBlockingQueue<MzmlBlock>(1);
			MzmlToDIAProducer producer=new MzmlToDIAProducer(unmarshaller, mzmlBlockQueue, parameters);
			
			// will be populated after we join back up. Since we're not looking
			// at it until after the join, we're safe to not have to worry about
			// concurrency.
			HashMap<Range, TFloatArrayList> retentionTimesByStripe=producer.getRetentionTimesByStripe();
			Thread[] threads;
			
			if (parameters.isDeconvoluteOverlappingWindows()) {
				BlockingQueue<MzmlBlock> deconvolutionBlockQueue=new ArrayBlockingQueue<MzmlBlock>(1);
				OverlapDeconvoluter deconvoluter=new OverlapDeconvoluter(parameters.getFragmentTolerance(), mzmlBlockQueue, deconvolutionBlockQueue);
				retentionTimesByStripe=deconvoluter.getRetentionTimesByStripe();
				MzmlToDIAConsumer consumer=new MzmlToDIAConsumer(deconvolutionBlockQueue, stripeFile);

				Logger.logLine("Converting "+mzMLFile.getName()+" ...");
				Thread producerThread=new Thread(producer);
				Thread deconvoluterThread=new Thread(deconvoluter);
				Thread consumerThread=new Thread(consumer);

				threads=new Thread[] {producerThread, deconvoluterThread, consumerThread};
				
			} else {
				MzmlToDIAConsumer consumer=new MzmlToDIAConsumer(mzmlBlockQueue, stripeFile);

				Logger.logLine("Converting "+mzMLFile.getName()+" ...");
				Thread producerThread=new Thread(producer);
				Thread consumerThread=new Thread(consumer);

				threads=new Thread[] {producerThread, consumerThread};
			}
			
			for (int i=0; i<threads.length; i++) {
				threads[i].start();
			}

			try {
				for (int i=0; i<threads.length; i++) {
					threads[i].join();
				}

				Logger.logLine("Finalizing "+diaFile.getName()+" ...");
				HashMap<Range, Float> dutyCycleMap=new HashMap<Range, Float>();
				for (Entry<Range, TFloatArrayList> entry : retentionTimesByStripe.entrySet()) {
					Range range=entry.getKey();
					TFloatArrayList rts=entry.getValue();
					float[] deltas=General.firstDerivative(rts.toArray());
					float averageDutyCycle=General.mean(deltas);
					dutyCycleMap.put(range, averageDutyCycle);
				}
				if (!parameters.isDDA()) {
					stripeFile.setRanges(dutyCycleMap);
				}

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
