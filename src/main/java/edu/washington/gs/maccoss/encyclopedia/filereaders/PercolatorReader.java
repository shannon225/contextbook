package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PercolatorReader {
	public static Pair<ArrayList<PercolatorPeptide>, Float> getPassingPeptidesFromTSV(File f, final float qValueThreshold, final boolean keepDecoys) {
		final ArrayList<PercolatorPeptide> data=new ArrayList<PercolatorPeptide>();
		final float[] pi0=new float[1];
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String psmID=row.get("PSMId");
				
				// PSMId is the first row, so any non-table data will get put into this column
				if (psmID.startsWith(PercolatorExecutor.PI_0_TAG)) {
					pi0[0]=Float.parseFloat(psmID.substring(PercolatorExecutor.PI_0_TAG.length()));
					return;
				}
				
				float qvalue=Float.parseFloat(row.get("q-value"));
				if (qvalue<qValueThreshold) {
					float posteriorErrorProb=Float.parseFloat(row.get("posterior_error_prob"));
					
					String proteinIds=row.get("proteinIds");
					if (keepDecoys||!proteinIds.startsWith(LibraryEntry.DECOY_STRING)) {
						data.add(new PercolatorPeptide(psmID, proteinIds, qvalue, posteriorErrorProb));
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
			Logger.errorLine("Percolator reading interrupted!");
			Logger.errorException(ie);
		}

		return new Pair<ArrayList<PercolatorPeptide>, Float>(data, pi0[0]);
	}
	
	/**
	 * THIS IS DANGEROUS BECAUSE THERE IS NO PARSIMONY APPLIED! CONSIDER JUST USING PARSIMONYPROTEINGROUPER
	 * @param f
	 * @param qValueThreshold
	 * @param keepDecoys
	 * @return
	 */
	public static ArrayList<PercolatorProteinGroup> getPassingProteinsFromTSV(File f, final float qValueThreshold, final boolean keepDecoys) {
		final ArrayList<PercolatorProteinGroup> data=new ArrayList<>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String proteinIDs=row.get("ProteinId");
				
				float qvalue=Float.parseFloat(row.get("q-value"));
				if (qvalue<qValueThreshold) {
					float posteriorErrorProb=Float.parseFloat(row.get("posterior_error_prob"));
					String peptideIDs=row.get("peptideIds");
					if (keepDecoys||!proteinIDs.startsWith(LibraryEntry.DECOY_STRING)) {
						data.add(new PercolatorProteinGroup(proteinIDs.split(PSMData.ACCESSION_TOKEN), peptideIDs.split(" "), qvalue, posteriorErrorProb));
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
			Logger.errorLine("Percolator reading interrupted!");
			Logger.errorException(ie);
		}

		return data;
	}
	
	public static ArrayList<ScoredObject<String>> getPassingPeptidesFromXML(File f, final float qValueThreshold) {
		ArrayList<ScoredObject<String>> data=readPercolator(f);
		
		ArrayList<ScoredObject<String>> thresholded=new ArrayList<ScoredObject<String>>();
		for (ScoredObject<String> scoredObject : data) {
			if (scoredObject.x<=qValueThreshold) {
				thresholded.add(scoredObject);
			}
		}
		return thresholded;
	}

	public static ArrayList<ScoredObject<String>> readPercolator(File f) {
		BufferedReader in=null;
		ArrayList<ScoredObject<String>> entryList=new ArrayList<ScoredObject<String>>();
		try {
			in=new BufferedReader(new FileReader(f));
			return readPercolator(in, f.getName());

		} catch (IOException ioe) {
			Logger.errorLine("I/O Error found reading Percolator file ["+f.getAbsolutePath()+"]");
			Logger.errorException(ioe);
			return entryList;
		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		}
	}
	
	public static ArrayList<ScoredObject<String>> readPercolator(String s, String fileName) {
		return readPercolator(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(s.getBytes(StandardCharsets.UTF_8)))), fileName);
	}
	
	public static ArrayList<ScoredObject<String>> readPercolator(InputStream s, String fileName) {
		return readPercolator(new BufferedReader(new InputStreamReader(s)), fileName);
	}
	
	public static ArrayList<ScoredObject<String>> readPercolator(BufferedReader in, String fileName) {
		ArrayList<ScoredObject<String>> entryList=new ArrayList<ScoredObject<String>>();
		try {
			boolean parsePeptides=false;
			float qValue=0.0f;
			
			String eachline;
			while ((eachline=in.readLine())!=null) {
				eachline=eachline.trim();
				if (eachline.length()==0) {
					continue;
				}
				if (parsePeptides) {
					if (eachline.startsWith("<q_value>")) {
						qValue=Float.parseFloat(eachline.substring(9, eachline.length()-10));
					} else if (eachline.startsWith("<psm_id>")) {
						String psmID=eachline.substring(8, eachline.length()-9);
						entryList.add(new ScoredObject<String>(qValue, psmID));
					} else if (eachline.startsWith("</peptides>")) {
						parsePeptides=false;
					}
				} else if (eachline.startsWith("<peptides>")) {
					parsePeptides=true;
				}
			}
			return entryList;

		} catch (IOException ioe) {
			Logger.errorLine("I/O Error found reading FASTA ["+fileName+"]");
			Logger.errorException(ioe);
			return entryList;
		} finally {
			if (in!=null) {
				try {
					in.close();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		}
	}
}
