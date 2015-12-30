package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

public class MzmlToDIAConverter {
	public static final String MZML_EXTENSION=".mzml";

	public static void main(String[] args) {
		Long time=System.currentTimeMillis();
		File xmlFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/20150708_Ecoli_0911_25x4mzDIA_500_600.mzML");
		File saveFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/20150708_Ecoli_0911_25x4mzDIA_500_600.dia");
		convert(xmlFile, saveFile);
		System.out.println((System.currentTimeMillis()-time)/1000f+" seconds");
	}

	public static StripeFileInterface getFile(File f) {
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
			return convert(f, diaFile);
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
			throw new EncyclopediaException("Error reading DIA file!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("Error reading DIA file!", sqle);
		}
	}

	static StripeFileInterface convert(File mzMLFile, File diaFile) {
		try {
			Logger.logLine("Indexing "+mzMLFile.getName()+" ...");
			StripeFile stripeFile=new StripeFile();
			stripeFile.openFile();

			MzMLUnmarshaller unmarshaller=new MzMLUnmarshaller(mzMLFile);
			stripeFile.setFileName(mzMLFile.getName(), unmarshaller.getMzMLId(), mzMLFile.getAbsolutePath());

			BlockingQueue<MzmlBlock> mzmlBlockQueue=new ArrayBlockingQueue<MzmlBlock>(1);
			MzmlToDIAProducer producer=new MzmlToDIAProducer(unmarshaller, mzmlBlockQueue);
			MzmlToDIAConsumer consumer=new MzmlToDIAConsumer(mzmlBlockQueue, stripeFile);

			Logger.logLine("Converting "+mzMLFile.getName()+" ...");
			Thread producerThread=new Thread(producer);
			Thread consumerThread=new Thread(consumer);
			producerThread.start();
			consumerThread.start();

			try {
				producerThread.join();
				consumerThread.join();

				Logger.logLine("Finalizing "+diaFile.getName()+" ...");
				HashMap<Range, Float> dutyCycleMap=new HashMap<Range, Float>();
				for (Entry<Range, TFloatArrayList> entry : producer.getRetentionTimesByStripe().entrySet()) {
					Range range=entry.getKey();
					TFloatArrayList rts=entry.getValue();
					float[] deltas=General.firstDerivative(rts.toArray());
					float averageDutyCycle=General.mean(deltas);
					dutyCycleMap.put(range, averageDutyCycle);
				}
				stripeFile.setRanges(dutyCycleMap);

				stripeFile.saveAsFile(diaFile);
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
			throw new EncyclopediaException("DIA writing SQL error!", sqle);
		}
	}
}
