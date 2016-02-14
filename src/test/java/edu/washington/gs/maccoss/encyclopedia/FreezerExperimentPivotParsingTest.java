package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class FreezerExperimentPivotParsingTest {
	
	public static void main(String[] args) throws Exception {
		File dir=new File("/Users/searleb/Documents/school/projects/freezer/quant/");
		HashSet<String> peptides=getValidPeptidesFromTSV(new File(dir, "consistent_peptides.txt"));
		System.out.println(peptides.size()+" total peptides");

		File[] files=dir.listFiles();
		long totalTime=0;
		int fileCount=0;
		
		int totalCount=0;
		for (File file : files) {
			if (file.getName().endsWith(".integration.txt")) {
				totalCount++;
			}
		}
		float[] timing=new float[] { 0.0f, 4.0f/24.0f, 12.0f/24.0f, 1f, 2f, 3f, 4f, 7f, 14f, 28f, 50f };

		TreeMap<String, float[]> table4c=new TreeMap<String, float[]>();
		TreeMap<String, float[]> tableRT=new TreeMap<String, float[]>();
		TreeMap<String, float[]> tableM20c=new TreeMap<String, float[]>();
		
		for (File file : files) {
			String filename=file.getName();
			if (filename.endsWith(".quant.txt")) {
				String timingToken=filename.substring("121115_bcs_hela_24mz_400_1000_".length()-1, filename.length()-".mzML.quant.txt".length());
				StringTokenizer st=new StringTokenizer(timingToken, "_");
				int count=st.countTokens();
				String type;
				int rep;
				float time;
				if (count==2) {
					// day 0 (0D_1)
					type="init";
					time=0.0f;
					st.nextToken();
					rep=Integer.parseInt(st.nextToken());
				} else if (count==3) {
					// non day 0 (4c_7D_2) (m20c_12H_1)
					type=st.nextToken();
					String timingString=st.nextToken();
					if (timingString.endsWith("D")) {
						time=Float.parseFloat(timingString.substring(0, timingString.length()-1));
					} else if (timingString.endsWith("H")) {
						time=Float.parseFloat(timingString.substring(0, timingString.length()-1))/24.0f;
					} else {
						throw new RuntimeException("Unknown file type for "+filename+": ["+timingToken+"]");
					}
					rep=Integer.parseInt(st.nextToken());
				} else {
					throw new RuntimeException("Unknown file type for "+filename+": ["+timingToken+"]");
				}
				
				int index=Arrays.binarySearch(timing, time);
				
				System.out.println(type+"\t"+rep+"\t"+time+"\t"+filename);
				ArrayList<ScoredObject<String>> data=getPeptidesFromTSV(file);
				
				TreeMap<String, float[]>[] dataGroups;
				if ("init".equals(type)) {
					dataGroups=new TreeMap[] {table4c, tableRT, tableM20c};
				} else if ("4c".equals(type)) {
					dataGroups=new TreeMap[] {table4c};
				} else if ("rt".equals(type)) {
					dataGroups=new TreeMap[] {tableRT};
				} else if ("m20c".equals(type)) {
					dataGroups=new TreeMap[] {tableM20c};
				} else {
					throw new RuntimeException("Unknown file type for "+filename+": ["+timingToken+"]");
				}
				
				for (TreeMap<String, float[]> map : dataGroups) {
					for (ScoredObject<String> peptide : data) {
						if (peptides.contains(peptide.y)) {
							float[] array=map.get(peptide.y);
							if (array==null) {
								array=new float[timing.length];
								map.put(peptide.y, array);
							}
							array[index]+=peptide.x;
						}
					}
				}
			}
		}
		
		writeTable(timing, new File(dir, "pivot_4c.txt"), table4c);
		writeTable(timing, new File(dir, "pivot_RT.txt"), tableRT);
		writeTable(timing, new File(dir, "pivot_m20c.txt"), tableM20c);
	}

	private static void writeTable(float[] timing, File file, TreeMap<String, float[]> map) throws FileNotFoundException, UnsupportedEncodingException {
		System.out.println("Writing "+file.getName()+" ("+map.size()+" peptides)");
		PrintWriter writer=new PrintWriter(file, "UTF-8");
		writer.print("peptide");
		for (float f : timing) {
			writer.print("\t"+f+"_Days");
		}
		writer.println();
		for (Entry<String, float[]> entry : map.entrySet()) {
			writer.print(entry.getKey());
			float[] array=entry.getValue();
			for (int i=0; i<array.length; i++) {
				writer.print("\t"+array[i]);
			}
			writer.println();
		}
		writer.close();
	}
	
	public static ArrayList<ScoredObject<String>> getPeptidesFromTSV(File f) {
		final ArrayList<ScoredObject<String>> data=new ArrayList<ScoredObject<String>>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String peptideModSeq=row.get("PeptideModSeq");
				float tic=Float.parseFloat(row.get("TIC"));
				data.add(new ScoredObject<String>(tic, peptideModSeq));
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
			Logger.errorLine("getPeptidesFromTSV reading interrupted!");
			Logger.errorException(ie);
		}

		return data;
	}

	public static HashSet<String> getValidPeptidesFromTSV(File f) {
		final HashSet<String> data=new HashSet<String>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				String peptideModSeq=row.get("PeptideModSeq");
				data.add(peptideModSeq);
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
			Logger.errorLine("getValidPeptidesFromTSV reading interrupted!");
			Logger.errorException(ie);
		}

		return data;
	}
}
