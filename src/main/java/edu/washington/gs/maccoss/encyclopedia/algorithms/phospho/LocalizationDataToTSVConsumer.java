package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector.OS;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

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
	
	public static HashMap<String, ModificationLocalizationData> readLocalizationFile(File f, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters parameters) {
		HashSet<String> passingPeptideModSeqs=new HashSet<>();
		for (PercolatorPeptide peptide : passingPeptides) {
			String peptideModSeq = PeptideUtils.getCorrectedMasses(peptide.getPeptideModSeq());
			passingPeptideModSeqs.add(peptideModSeq);
		}

		final PeptideModification modification=parameters.getLocalizingModification().get();
		final HashMap<String, ModificationLocalizationData> result=new HashMap<String, ModificationLocalizationData>();
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String peptideModSeq=row.get("peptideModSeq");
				if (passingPeptideModSeqs.contains(peptideModSeq)) {
					float localizationScore=Float.parseFloat(row.get("localizationScore"));
					boolean isSiteSpecific=Boolean.parseBoolean(row.get("isSiteSpecific"));
					boolean isLocalized=Boolean.parseBoolean(row.get("isLocalized"));
					ModificationLocalizationData prev=result.get(peptideModSeq);
					
					boolean moreSiteSpecific = prev==null||(isSiteSpecific&&!prev.isSiteSpecific());
					boolean notLessSiteSpecific = prev!=null&&(isSiteSpecific||(!prev.isSiteSpecific()));
					boolean higherScoring = prev==null||(notLessSiteSpecific&&prev.getLocalizationScore()<localizationScore);
					
					if (moreSiteSpecific||higherScoring) {
						try {
							AmbiguousPeptideModSeq localizationPeptideModSeq=AmbiguousPeptideModSeq.getAmbiguousPeptideModSeq(peptideModSeq, modification);
							float retentionTimeApexInSeconds=Float.parseFloat(row.get("retentionTimeApexInSeconds"));
							int numberOfMods=Integer.parseInt(row.get("numberOfMods"));
							FragmentIon[] localizingIons;
							try {
								localizingIons=FragmentIon.fromArchiveString(row.get("localizingIons"));
							} catch (Exception e) {
								Logger.errorLine("Error parsing localization ions for "+peptideModSeq+"from ["+row.get("localizingIons")+"], skipping ions but keeping peptide.");
								localizingIons=new FragmentIon[0];
							}
							float localizingIntensity=Float.parseFloat(row.get("localizingIntensity"));
							float totalIntensity=Float.parseFloat(row.get("totalIntensity"));
							
							ModificationLocalizationData data=new ModificationLocalizationData(localizationPeptideModSeq, retentionTimeApexInSeconds, localizationScore, numberOfMods, isSiteSpecific, isLocalized, localizingIons, localizingIntensity, totalIntensity);
							result.put(peptideModSeq, data);
						} catch (Exception e) {
							Logger.errorLine("Error parsing localization data for "+peptideModSeq+", skipping this peptide! ("+e.getMessage()+")");
							e.printStackTrace();
						}
					}
				}
			}
		};
		
		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, f, "\t", 1);
		TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
		} catch (InterruptedException ie) {
			Logger.errorLine("Localization reading interrupted!");
			Logger.errorException(ie);
		}
		return result;
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
					writer.print("peptideModSeq\tlocalizationPeptideModSeq\tretentionTimeApexInSeconds\tlocalizationScore\tnumberOfMods\tisSiteSpecific\tisLocalized\tlocalizingIons\tlocalizingIntensity\ttotalIntensity");
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
				writer.print("\t"+result.isLocalized());
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
