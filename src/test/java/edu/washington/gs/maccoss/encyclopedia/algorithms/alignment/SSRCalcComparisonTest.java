package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SSRCalc;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class SSRCalcComparisonTest {
	public static void main(String[] args) {
		File actualFile = new File("/Users/searleb/Documents/oleg/from_me/peptides_for_Oleg.csv");
		File ssrFile = new File("/Users/searleb/Documents/oleg/from_me/BRS_PRED_HIDB.txt");
		File prositFile = new File("/Users/searleb/Documents/oleg/from_me/peptide_sequences_prosit.csv");
		File resnetFile = new File("/Users/searleb/Documents/oleg/from_me/peptide_sequences_resnet.csv");
		File oldssrFile = new File("/Users/searleb/Documents/oleg/from_me/oldSSRCalc");
		
		TObjectFloatHashMap<String> actual=getData(actualFile, "PeptideModSeq", "RTInSeconds", ",", 1/60f);
		TObjectFloatHashMap<String> ssrcalc=getData(ssrFile, "Seq", "HI2018", "\t", 1);
		TObjectFloatHashMap<String> prositcalc=getData(prositFile, "PeptideModSeq", "RTInSeconds", ",", 1/60f);
		TObjectFloatHashMap<String> resnetcalc=getData(resnetFile, "PeptideSeq", "RTInSeconds", ",", 1/60f);
		
		ArrayList<XYPoint> actualSSRCalcPair=new ArrayList<>();
		ArrayList<XYPoint> actualPrositPair=new ArrayList<>();
		ArrayList<XYPoint> actualResnetPair=new ArrayList<>();
		ArrayList<XYPoint> actualOldSSRCalcPair=new ArrayList<>();
		actual.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float actual) {
				float ssr=ssrcalc.get(pep);
				float prosit=prositcalc.get(pep);
				float resnet=resnetcalc.get(pep);
				
				if (ssr!=0.0f&&prosit!=0.0f&&resnet!=0.0f) {
					actualSSRCalcPair.add(new RTRTPoint(ssr, actual, false, pep));
					actualPrositPair.add(new RTRTPoint(prosit, actual, false, pep));
					actualResnetPair.add(new RTRTPoint(resnet, actual, false, pep));
					actualOldSSRCalcPair.add(new RTRTPoint(SSRCalc.getHydrophobicity(pep), actual, false, pep));
					//System.out.println(resnet+"\t"+SSRCalc.getHydrophobicity(pep));
				}
				return true;
			}
		});
		
		RetentionTimeFilter filter=RetentionTimeFilter.getFilter(actualSSRCalcPair);
		filter.plot(actualSSRCalcPair, Optional.ofNullable(ssrFile));

		RetentionTimeFilter prositFilter=RetentionTimeFilter.getFilter(actualPrositPair);
		prositFilter.plot(actualPrositPair, Optional.ofNullable(prositFile));

		RetentionTimeFilter resnetFilter=RetentionTimeFilter.getFilter(actualResnetPair);
		resnetFilter.plot(actualResnetPair, Optional.ofNullable(resnetFile));

		RetentionTimeFilter oldSSRCalcFilter=RetentionTimeFilter.getFilter(actualOldSSRCalcPair);
		oldSSRCalcFilter.plot(actualOldSSRCalcPair, Optional.ofNullable(oldssrFile));
	}

	public static TObjectFloatHashMap<String> getData(File f, String pepHeader, String rtHeader, String delim, float multiplier) {
		TObjectFloatHashMap<String> map=new TObjectFloatHashMap<>();
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String pep=row.get(pepHeader);
				pep=pep.replace("[+57.021464]", "");
				
				String s=row.get(rtHeader);
				float predicted=Float.parseFloat(s)*multiplier;
				map.put(pep, predicted);
			}
			
			@Override
			public void cleanup() {
			}
		};

		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, f, delim, 1);
		TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
		} catch (InterruptedException ie) {
			Logger.errorLine("Error reading data!");
			Logger.errorException(ie);
		}

		return map;
	}
}
