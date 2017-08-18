package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
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
import gnu.trove.map.hash.TObjectFloatHashMap;

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
	
	public static class AmbiguouslyModifiedPeptide implements Comparable<AmbiguouslyModifiedPeptide> {
		private final String peptideSeq;
		private final double precursorMz;
		private final byte precursorCharge;
		public AmbiguouslyModifiedPeptide(String peptideSeq, double precursorMz, byte precursorCharge) {
			this.peptideSeq=peptideSeq;
			this.precursorMz=precursorMz;
			this.precursorCharge=precursorCharge;
		}
		public AmbiguouslyModifiedPeptide(PercolatorPeptide peptide, AminoAcidConstants aaConstants) {
			String peptideModSeq=peptide.getPeptideModSeq();
			peptideSeq=PeptideUtils.getPeptideSeq(peptideModSeq);
			precursorCharge=peptide.getPrecursorCharge();
			precursorMz=aaConstants.getChargedMass(peptideModSeq, precursorCharge);
		}
		@Override
		public int compareTo(AmbiguouslyModifiedPeptide o) {
			if (o==null) return 1;
			int c=peptideSeq.compareTo(o.peptideSeq);
			if (c!=0) return c;
			c=Double.compare(precursorMz, o.precursorMz);
			if (c!=0) return c;
			c=Byte.compare(precursorCharge, o.precursorCharge);
			return c;
		}
		
		@Override
		public boolean equals(Object obj) {
			return compareTo((AmbiguouslyModifiedPeptide)obj)==0;
		}
		@Override
		public int hashCode() {
			return peptideSeq.hashCode()+Double.hashCode(precursorMz)+precursorCharge;
		}
	}
	
	public static HashMap<String, ModificationLocalizationData> readLocalizationFile(File f, ArrayList<PercolatorPeptide> passingPeptides, SearchParameters parameters) {
		ArrayList<PercolatorPeptide> cloned=new ArrayList<PercolatorPeptide>(passingPeptides);
		Collections.sort(cloned, PercolatorPeptide.scoreComparator);
		
		HashMap<String, PercolatorPeptide> passingPeptideModSeqs=new HashMap<>();
		for (PercolatorPeptide peptide : cloned) {
			passingPeptideModSeqs.put(peptide.getPeptideModSeq(), peptide);
		}
		final AminoAcidConstants aaConstants=parameters.getAAConstants();
		final PeptideModification modification=parameters.getLocalizingModification().get();
		final HashMap<String, ModificationLocalizationData> result=new HashMap<String, ModificationLocalizationData>();
		final HashSet<AmbiguouslyModifiedPeptide> previouslyDetected=new HashSet<>();
		final TObjectFloatHashMap<AmbiguouslyModifiedPeptide> notAnnotated=new TObjectFloatHashMap<>();


		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String peptideModSeq=row.get("peptideModSeq");
				float totalIntensity=Float.parseFloat(row.get("totalIntensity"));
				AmbiguouslyModifiedPeptide ambigous=new AmbiguouslyModifiedPeptide(passingPeptideModSeqs.remove(peptideModSeq), aaConstants);
				
				if (passingPeptideModSeqs.containsKey(peptideModSeq)) {
					previouslyDetected.add(ambigous);
					float localizationScore=Float.parseFloat(row.get("localizationScore"));
					boolean isSiteSpecific=Boolean.parseBoolean(row.get("isSiteSpecific"));
					boolean isLocalized=Boolean.parseBoolean(row.get("isLocalized"));
					ModificationLocalizationData prev=result.get(peptideModSeq);
					String ambiPeptideModSeq=row.get("localizationPeptideModSeq");
					AmbiguousPeptideModSeq localizationPeptideModSeq=AmbiguousPeptideModSeq.getAmbiguousPeptideModSeq(ambiPeptideModSeq, modification);
					
					boolean moreSiteSpecific = prev==null||(isSiteSpecific&&!prev.isSiteSpecific());
					boolean notLessSiteSpecific = prev==null||localizationPeptideModSeq.getAmbiguityValue()>=prev.getLocalizationPeptideModSeq().getAmbiguityValue();
					boolean higherScoring = prev==null||(notLessSiteSpecific&&prev.getLocalizationScore()<localizationScore);
					boolean isCompletelyAmbiguous=AmbiguousPeptideModSeq.isCompletelyAmbiguous(localizationPeptideModSeq, modification);
					
					if (moreSiteSpecific||higherScoring) {
						try {
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
							
							ModificationLocalizationData data=new ModificationLocalizationData(localizationPeptideModSeq, retentionTimeApexInSeconds, localizationScore, numberOfMods, isSiteSpecific, isLocalized, isCompletelyAmbiguous, localizingIons, localizingIntensity, totalIntensity);
							result.put(peptideModSeq, data);
						} catch (Exception e) {
							Logger.errorLine("Error parsing localization data for "+peptideModSeq+", skipping this peptide! ("+e.getMessage()+")");
							e.printStackTrace();
						}
					}
				} else {
					// NOTE: relies on empty values being 0
					notAnnotated.put(ambigous, Math.max(notAnnotated.get(ambigous), totalIntensity));
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

		for (PercolatorPeptide peptide : cloned) {
			AmbiguouslyModifiedPeptide ambigPeptide=new AmbiguouslyModifiedPeptide(peptide, parameters.getAAConstants());
			float totalIntensity=notAnnotated.get(ambigPeptide);

			if (!previouslyDetected.contains(ambigPeptide)) {
				previouslyDetected.add(ambigPeptide);
				
				AmbiguousPeptideModSeq localizationPeptideModSeq=AmbiguousPeptideModSeq.getFullyAmbiguous(peptide.getPeptideModSeq(), modification, aaConstants);
				
				ModificationLocalizationData data=new ModificationLocalizationData(localizationPeptideModSeq, peptide.getRT(), 0.0f, localizationPeptideModSeq.getNumModifications(), false, false, true, new FragmentIon[0], 0.0f, totalIntensity);
				result.put(peptide.getPeptideModSeq(), data);
			}
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
