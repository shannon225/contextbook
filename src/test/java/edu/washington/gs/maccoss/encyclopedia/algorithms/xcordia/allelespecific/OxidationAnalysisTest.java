package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.XCorDIA;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeFilter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideDatabase;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class OxidationAnalysisTest {
	public static void main(String[] args) throws Exception {
		HashMap<String, String> defaults=PecanParameterParser.getDefaultParameters();
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(defaults);

		File diaFile=new File("/Users/searleb/Documents/school/xcordia/hela_prm/2017dec27_normal_dia_6b_rep1.dia");
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		ArrayList<Range> ranges=new ArrayList<>(stripefile.getRanges().keySet());
		Collections.sort(ranges);
		
		PeptideDatabase targets=new PeptideDatabase();
		PeptideDatabase possibleTargets=new PeptideDatabase();
		ArrayList<PeptideData> widePeptides=readData(new File("/Users/searleb/Documents/school/xcordia/hela_prm/unique_peptides.txt"));
		ArrayList<PeptideData> narrowPeptides=readData(new File("/Users/searleb/Documents/school/xcordia/hela_prm/librarypeptides.txt"));
		

		TObjectFloatHashMap<String> widePeptideRTs=new TObjectFloatHashMap<>();
		for (PeptideData peptideData : widePeptides) {
			targets.add(peptideData);
			possibleTargets.add(peptideData);
			widePeptideRTs.put(peptideData.getSequence(), peptideData.rt);
		}
		
		TObjectFloatHashMap<String> libraryPeptideRTs=new TObjectFloatHashMap<>();
		for (PeptideData peptideData : narrowPeptides) {
			possibleTargets.add(peptideData);
			libraryPeptideRTs.put(peptideData.getSequence(), peptideData.rt);
		}

		ArrayList<XYPoint> rtPoints=new ArrayList<>();
		for (Range range : ranges) {
			HashSet<FastaPeptideEntry> peptides=XCorDIA.getPeptidesInRange(parameters, targets, range);
			SimilarPeptideBinner binner=new SimilarPeptideBinner();
			ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(peptides);
			
			for (ArrayList<FastaPeptideEntry> list : bins) {
				if (list.size()==1) {
					PeptideData peptideData=(PeptideData)list.get(0);
					float wideRT = widePeptideRTs.get(peptideData.getSequence());
					float narrowRT = libraryPeptideRTs.get(peptideData.getSequence());
					if (wideRT>20&&narrowRT>20) {
						rtPoints.add(new XYPoint(narrowRT, wideRT));
					}
				}
			}
		}
		RetentionTimeFilter filter=RetentionTimeFilter.getFilter(rtPoints, "Library", "Wide");
		filter.plot(rtPoints, Optional.empty());
		
		int betterAltDeltas=0;
		for (Range range : ranges) {
			HashSet<FastaPeptideEntry> peptides=XCorDIA.getPeptidesInRange(parameters, possibleTargets, range);
			SimilarPeptideBinner binner=new SimilarPeptideBinner();
			ArrayList<ArrayList<FastaPeptideEntry>> bins=binner.binPeptides(peptides);
			
			for (ArrayList<FastaPeptideEntry> list : bins) {
				
				if (list.size()>1) {
					for (FastaPeptideEntry entry : list) {
						float thisRT=widePeptideRTs.get(entry.getSequence());
						if (thisRT==0) continue; // this peptide wasn't detected in the wide search
						
						thisRT=filter.getXValue(thisRT); // map to Library
						float libRT=libraryPeptideRTs.get(entry.getSequence());
						if (libRT==0.0f) libRT=-1000.0f; 
						float deltaRT=libRT-thisRT;
	
						FastaPeptideEntry bestAlt=null;
						float altDeltaRT=+1000.0f;
						float altLibRT=0.0f;
						for (FastaPeptideEntry alt : list) {
							if (alt!=entry) {
								float thisAltLibRT=libraryPeptideRTs.get(alt.getSequence());
								if (thisAltLibRT==0.0f) thisAltLibRT=-1000.0f; 
								float thisDeltaRT=thisAltLibRT-thisRT;
								if (Math.abs(thisDeltaRT)<Math.abs(altDeltaRT)) {
									altDeltaRT=thisDeltaRT;
									bestAlt=alt;
									altLibRT=thisAltLibRT;
								}
							}
						}
						String altSequence = bestAlt==null?"none":bestAlt.getSequence();
						if (bestAlt!=null&&list.size()>1&&Math.abs(altDeltaRT)<Math.abs(deltaRT)) {
							System.out.println(range+"\t"+entry.getSequence()+"\t"+thisRT+"\t"+list.size()+"\t"+deltaRT+"\t"+altDeltaRT+"\t"+altSequence+"\t"+altLibRT);
							betterAltDeltas++;
						}
					}
				}
			}
		}
		System.out.println("total peptides:"+targets.size()+"\tbetter RTs with alts:"+betterAltDeltas);
	}
	
	private static ArrayList<PeptideData> readData(File f) {
		ArrayList<PeptideData> peptides=new ArrayList<>();
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String sequence=row.get("peptidemodseq");
				String filename=row.get("sourcefile");
				float rt=Float.parseFloat(row.get("rtinmin"));
				peptides.add(new PeptideData(filename, sequence, rt));
			}
			
			@Override
			public void cleanup() {
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
			Logger.errorLine("Percolator reading interrupted!");
			Logger.errorException(ie);
		}

		return peptides;
	}

	private static class PeptideData extends FastaPeptideEntry {
		private final float rt;
		public PeptideData(String filename, String sequence, float rt) {
			super(filename, "", sequence);
			this.rt=rt;
		}
	}
}
