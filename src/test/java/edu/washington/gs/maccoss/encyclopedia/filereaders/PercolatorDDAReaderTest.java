package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class PercolatorDDAReaderTest {

	private static final String DECOY = "DECOY_";
	private static final String TARGET = "TARGET_";

	public static void main(String[] args) {
		File ft=new File("/Users/searleb/Documents/damien/dda_library_search/msfragger_hela/perc_results/23aug2017_hela_serum_timecourse_pool_dda_002_percolator_target_psms.tsv");
		File fd=new File("/Users/searleb/Documents/damien/dda_library_search/msfragger_hela/perc_results/23aug2017_hela_serum_timecourse_pool_dda_002_percolator_decoy_psms.tsv");
		AtomicInteger rowNumber=new AtomicInteger(0);
		
		HashSet<String> observedPeptides=new HashSet<>();
		ArrayList<ScoredObject<String>> scoredPeptides=new ArrayList<>();
		TableParserMuscle targetMuscle = getMuscle(rowNumber, true, observedPeptides, scoredPeptides);
		TableParserMuscle decoyMuscle = getMuscle(rowNumber, false, observedPeptides, scoredPeptides);

		readFile(ft, targetMuscle);
		System.out.println(rowNumber+" vs "+scoredPeptides.size());
		readFile(fd, decoyMuscle);
		System.out.println(rowNumber+" vs "+scoredPeptides.size());
		
		Collections.sort(scoredPeptides);
		Collections.reverse(scoredPeptides);
		
		int totalTargets=0;
		int totalDecoys=0;
		for (ScoredObject<String> pep : scoredPeptides) {
			if (pep.y.startsWith(TARGET)) totalTargets++;
			if (pep.y.startsWith(DECOY)) totalDecoys++;
			float fdr = totalDecoys/(float)totalTargets;
			System.out.println(totalTargets+"\t"+fdr);
			if (fdr>=0.01f) break;
		}
	}

	private static TableParserMuscle getMuscle(AtomicInteger rowNumber, final boolean isTarget,
			HashSet<String> observedPeptides, ArrayList<ScoredObject<String>> scoredPeptides) {
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				rowNumber.incrementAndGet();
				float score=Float.parseFloat(row.get("score"));
				String peptide=row.get("peptide");
				peptide=peptide.subSequence(peptide.indexOf('.')+1, peptide.lastIndexOf('.')).toString();
				peptide=PeptideUtils.getPeptideSeq(peptide);
				
				if (observedPeptides.contains(peptide)) return;
				observedPeptides.add(peptide);
				String id=(isTarget?TARGET:DECOY)+peptide;
						
				scoredPeptides.add(new ScoredObject<String>(score, id));
			}
			
			@Override
			public void cleanup() {
			}
		};
		return muscle;
	}

	private static void readFile(File ft, TableParserMuscle muscle) {
		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, ft, "\t", 1);
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
	}
}
