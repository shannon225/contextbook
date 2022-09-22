package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class LibraryRTComparisonTest {
	public static void main(String[] args) throws Exception {
		File cidPhosphoFile = new File("/Users/searleb/Documents/phospho/oleg_alignment/CID_phospho_rts.csv");
		File hcdPhosphoFile = new File("/Users/searleb/Documents/phospho/oleg_alignment/HCD_phospho_rts.csv");
		File hi2018File = new File("/Users/searleb/Documents/phospho/oleg_alignment/Oleg_10KPHOS_reformatted_DBW211015.txt");
		File reportFile = new File("/Users/searleb/Documents/phospho/oleg_alignment/phospho_combined_libraries.tsv");
		
		TObjectFloatHashMap<String> cidPhospho=getData(cidPhosphoFile, "PeptideModSeq", "rtinmin", ",", 1);
		TObjectFloatHashMap<String> hcdPhospho=getData(hcdPhosphoFile, "PeptideModSeq", "rtinmin", ",", 1);
		TObjectFloatHashMap<String> hi2018=getData(hi2018File, "peptide", "hydrophobicity", "\t", 1);
		
		ArrayList<XYPoint> cidPhosphoPair=new ArrayList<>();
		ArrayList<XYPoint>hcdPhosphoPair=new ArrayList<>();
		hi2018.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float h) {
				float b=cidPhospho.get(pep);
				float p=hcdPhospho.get(pep);
				
				if (h!=0.0f&&b!=0.0f) {
					cidPhosphoPair.add(new PeptideXYPoint(b, h, false, pep));
				}
				if (h!=0.0f&&p!=0.0f) {
					hcdPhosphoPair.add(new PeptideXYPoint(p, h, false, pep));
				}
				return true;
			}
		});
		
		FileWriter writer=new FileWriter(reportFile);
		final PrintWriter out=new PrintWriter(writer);
		out.println("PeptideModSeq\tHydrophobicity\tType");

		hi2018.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float h) {
				out.println(pep+"\t"+h+"\t"+"h18");
				return true;
			}
		});
		
		final RetentionTimeFilter cidPhosphoFilter=RetentionTimeFilter.getFilter(cidPhosphoPair);
		cidPhosphoFilter.plot(cidPhosphoPair, Optional.ofNullable(cidPhosphoFile));

		cidPhospho.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float b) {
				out.println(pep+"\t"+cidPhosphoFilter.getYValue(b)+"\t"+"cid");
				return true;
			}
		});

		final RetentionTimeFilter hcdPhosphoFilter=RetentionTimeFilter.getFilter(hcdPhosphoPair);
		hcdPhosphoFilter.plot(hcdPhosphoPair, Optional.ofNullable(hcdPhosphoFile));

		hcdPhospho.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float p) {
				out.println(pep+"\t"+hcdPhosphoFilter.getYValue(p)+"\t"+"hcd");
				return true;
			}
		});
		
		out.close();
	}
	public static void main2(String[] args) throws Exception {
		File bruderer2017File = new File("/Users/searleb/Documents/oleg/combined_libs/Bruderer2017_median-iRT_DBW211005.tsv");
		File proteomeToolsFile = new File("/Users/searleb/Documents/oleg/combined_libs/proteometools_rts_with_ox.csv");
		File hi2018File = new File("/Users/searleb/Documents/oleg/combined_libs/WHI2018-BRIAN400K.tsv");
		File reportFile = new File("/Users/searleb/Documents/oleg/combined_libs/combined_libraries.tsv");
		
		TObjectFloatHashMap<String> bruderer2017=getData(bruderer2017File, "PeptideModSeq", "iRT", "\t", 1);
		TObjectFloatHashMap<String> proteomeTools=getData(proteomeToolsFile, "PeptideModSeq", "iRT", ",", 1);
		TObjectFloatHashMap<String> hi2018=getData(hi2018File, "peptide", "hydrophobicity", "\t", 1);
		
		ArrayList<XYPoint> hiVSbruderer2017Pair=new ArrayList<>();
		ArrayList<XYPoint> hiVSproteomeToolsPair=new ArrayList<>();
		hi2018.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float h) {
				float b=bruderer2017.get(pep);
				float p=proteomeTools.get(pep);
				
				if (h!=0.0f&&b!=0.0f) {
					hiVSbruderer2017Pair.add(new PeptideXYPoint(b, h, false, pep));
				}
				if (h!=0.0f&&p!=0.0f) {
					hiVSproteomeToolsPair.add(new PeptideXYPoint(p, h, false, pep));
				}
				return true;
			}
		});
		
		FileWriter writer=new FileWriter(reportFile);
		final PrintWriter out=new PrintWriter(writer);
		out.println("PeptideModSeq\tHydrophobicity\thbp");

		hi2018.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float h) {
				out.println(pep+"\t"+h+"\t"+"h");
				return true;
			}
		});
		
		final RetentionTimeFilter bruderer2017Filter=RetentionTimeFilter.getFilter(hiVSbruderer2017Pair);
		bruderer2017Filter.plot(hiVSbruderer2017Pair, Optional.ofNullable(bruderer2017File));

		bruderer2017.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float b) {
				out.println(pep+"\t"+bruderer2017Filter.getYValue(b)+"\t"+"b");
				return true;
			}
		});

		final RetentionTimeFilter proteomeToolsFilter=RetentionTimeFilter.getFilter(hiVSproteomeToolsPair);
		proteomeToolsFilter.plot(hiVSproteomeToolsPair, Optional.ofNullable(proteomeToolsFile));

		proteomeTools.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String pep, float p) {
				out.println(pep+"\t"+proteomeToolsFilter.getYValue(p)+"\t"+"p");
				return true;
			}
		});
		
		out.close();
	}

	public static TObjectFloatHashMap<String> getData(File f, String pepHeader, String rtHeader, String delim, float multiplier) {
		TObjectFloatHashMap<String> map=new TObjectFloatHashMap<>();
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String pep=row.get(pepHeader);
				pep=pep.replace("[+57.021464]", "");
				pep=pep.replace("[79.966331]", "[+79.966331]");
				
				
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
