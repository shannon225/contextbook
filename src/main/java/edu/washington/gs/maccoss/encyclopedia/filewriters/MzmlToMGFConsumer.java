package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MSMSBlock;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlSAXToMSMSProducer;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;

public class MzmlToMGFConsumer implements Runnable {
	private final BlockingQueue<MSMSBlock> mzmlBlockQueue;
	private final File f;

	public MzmlToMGFConsumer(BlockingQueue<MSMSBlock> mzmlBlockQueue, File f) {
		this.mzmlBlockQueue=mzmlBlockQueue;
		this.f=f;
	}

	@Override
	public void run() {
		PrintWriter writer=null;
		try {
			writer=new PrintWriter(f, "UTF-8");
			while (true) {
				MSMSBlock block=mzmlBlockQueue.take();
				if (MSMSBlock.POISON_BLOCK==block) break;
				
				for (FragmentScan stripe : block.getFragmentScans()) {
					byte charge=stripe.getCharge();
					if (charge==0) continue;
					
					writer.println("BEGIN IONS");
					if (charge>0) {
						writer.print("PEPMASS=");
						writer.println(MassConstants.getPeptideMass(stripe.getPrecursorMZ(), charge));
						
						writer.print("CHARGE=");
						writer.print(charge);
						writer.println("+");
					}
					writer.print("TITLE=");
					writer.println(stripe.getSpectrumName());
					
					double[] masses=stripe.getMassArray();
					float[] intensities=stripe.getIntensityArray();
					for (int i=0; i<masses.length; i++) {
						writer.print(masses[i]);
						writer.print('\t');
						writer.println(intensities[i]);
					}

					writer.println("END IONS");
				}
			}
			

		} catch (InterruptedException e) {
			Logger.logException(e);
		} catch (FileNotFoundException e) {
			Logger.logException(e);
		} catch (UnsupportedEncodingException e) {
			Logger.logException(e);
		} finally {
			if (writer!=null) {
				writer.close();
			}
		}
	}
	
	public static void main(String[] args) {
		File mzMLFile=new File("/Users/searleb/Documents/damien/dda_library_search/hela/23aug2017_hela_serum_timecourse_pool_dda_003.mzML");
		File mgfFile=new File("/Users/searleb/Documents/damien/dda_library_search/hela/23aug2017_hela_serum_timecourse_pool_dda_003.mgf");
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		convertSAX(mzMLFile, mgfFile, parameters);
	}

	static void convertSAX(File mzMLFile, File mgfFile, SearchParameters parameters) {
		Logger.logLine("Indexing "+mzMLFile.getName()+" ...");

		BlockingQueue<MSMSBlock> mzmlBlockQueue=new ArrayBlockingQueue<MSMSBlock>(1);
		MzmlSAXToMSMSProducer producer=new MzmlSAXToMSMSProducer(mzMLFile, 0, mzmlBlockQueue, parameters);

		Thread[] threads;
		MzmlToMGFConsumer consumer=new MzmlToMGFConsumer(mzmlBlockQueue, mgfFile);

		Logger.logLine("Converting "+mzMLFile.getName()+" ...");
		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);

		threads=new Thread[] { producerThread, consumerThread };

		for (int i=0; i<threads.length; i++) {
			threads[i].start();
		}

		try {
			for (int i=0; i<threads.length; i++) {
				threads[i].join();
			}
			Logger.logLine("Finished writing "+mgfFile.getName()+"!");

		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		}

	}
}
