package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;

public class LocalizationDataToTSVConsumer implements Runnable {
	private final OS os=OSDetector.getOS();
	
	private final BlockingQueue<ModificationLocalizationData> resultsQueue;
	private final PrintWriter writer;
	private volatile int numberProcessed=0;

	public LocalizationDataToTSVConsumer(File outputFile, BlockingQueue<ModificationLocalizationData> resultsQueue) {
		this.resultsQueue=resultsQueue;
		try {
			writer=new PrintWriter(outputFile, "UTF-8");
			System.out.println("Constructing writer for "+outputFile.getAbsolutePath());
		} catch (FileNotFoundException e) {
			throw new EncyclopediaException("Error setting up output file: "+outputFile.getAbsolutePath(), e);
		} catch (UnsupportedEncodingException e) {
			throw new EncyclopediaException("Error setting up output file: "+outputFile.getAbsolutePath(), e);
		}
	}

	public void close() {
		writer.flush();
		writer.close();
	}

	public int getNumberProcessed() {
		return numberProcessed;
	}
	
	@Override
	public void run() {
		boolean printedHeader=false; 
		try {
			while (true) {
				ModificationLocalizationData result=resultsQueue.take();
				
				if (ModificationLocalizationData.POISON_RESULT==result) break;
				numberProcessed++;
				
				if (!printedHeader) {
					/**
						private final AmbiguousPeptideModSeq localizationPeptideModSeq;
						private final float retentionTimeApexInSeconds;
						private final float localizationScore;
						private final int numberOfMods;
						private final boolean isSiteSpecific;
						private final FragmentIon[] localizingIons;
						private final float localizingIntensity;
						private final float totalIntensity;
					 */
					writer.print("peptideModSeq\tlocalizationPeptideModSeq\tretentionTimeApexInSeconds\tlocalizationScore\tnumberOfMods\tisSiteSpecific\tlocalizingIons\tlocalizingIntensity\ttotalIntensity");
					// Percolator assumes linux line endings on Mac!
					switch (os) {
						case MAC:
							writer.print("\n");
							break;

						default:
							writer.println();
							break;
					}
					printedHeader=true;
				}
				writer.print(result.getLocalizationPeptideModSeq().getPeptideModSeq());
				writer.print("\t"+result.getLocalizationPeptideModSeq().getPeptideAnnotation());
				writer.print("\t"+result.getRetentionTimeApexInSeconds());
				writer.print("\t"+result.getLocalizationScore());
				writer.print("\t"+result.getNumberOfMods());
				writer.print("\t"+result.isSiteSpecific());
				writer.print("\t"+FragmentIon.toArchiveString(result.getLocalizingIons()));
				writer.print("\t"+result.getLocalizingIntensity());
				writer.print("\t"+result.getTotalIntensity());
				
				// Percolator assumes linux line endings on Mac!
				switch (os) {
					case MAC:
						writer.print("\n");
						break;

					default:
						writer.println();
						break;
				}
			}
		} catch (InterruptedException ie) {
			Logger.errorLine("DIA writing interrupted!");
			Logger.errorException(ie);
		}
	}
}
