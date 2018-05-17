package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import org.apache.commons.math3.distribution.HypergeometricDistribution;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import gnu.trove.map.hash.TObjectIntHashMap;

public class FishersExactTest {
	public static void main(String[] args) {
		File[] files=new File("/Users/searleb/Downloads/emma/").listFiles();
		for (File f : files) {

			ArrayList<String> groups=new ArrayList<>();
			float[] groupTotals=new float[10];
			TObjectIntHashMap<String> headerMap=new TObjectIntHashMap<>();

			TableParserMuscle muscle=new TableParserMuscle() {
				@Override
				public void processRow(Map<String, String> row) {
					if (groups.size()==0) {
						for (String header : row.keySet()) {
							if ("protid".equals(header)) continue;

							String group=header.substring(0, 2);
							int groupIndex=-1;
							for (int i=0; i<groups.size(); i++) {
								if (groups.get(i).equals(group)) {
									groupIndex=i;
									break;
								}
							}
							if (groupIndex==-1) {
								groupIndex=groups.size();
								groups.add(group);
							}
							headerMap.put(header, groupIndex);
						}
						
					}
					for (Entry<String, String> entry : row.entrySet()) {
						if ("protid".equals(entry.getKey())) continue;
						
						int index=headerMap.get(entry.getKey());
						groupTotals[index]+=Float.parseFloat(entry.getValue());
					}

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

			System.out.println(f.getName()+": Comparing between "+General.toString(groups.toArray())+"\tpvalue");

			muscle=new TableParserMuscle() {
				@Override
				public void processRow(Map<String, String> row) {
					String protein=null;
					float[] data=new float[groups.size()];
					
					for (Entry<String, String> entry : row.entrySet()) {
						if ("protid".equals(entry.getKey())) {;
							protein=entry.getValue();
							continue;
						}
						
						int index=headerMap.get(entry.getKey());
						
						data[index]+=Float.parseFloat(entry.getValue());
					}
					
					double pvalue=getFishersExactPvalue(groupTotals[0], groupTotals[1], data[0], data[1]);
					System.out.println(protein+"\t"+pvalue);
				}
				
				@Override
				public void cleanup() {
				}
			};

			blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
			producer=new TableParserProducer(blockingQueue, f, "\t", 1);
			consumer=new TableParserConsumer(blockingQueue, muscle);

			producerThread=new Thread(producer);
			consumerThread=new Thread(consumer);
			producerThread.start();
			consumerThread.start();

			try {
				producerThread.join();
				consumerThread.join();
			} catch (InterruptedException ie) {
				Logger.errorLine("Percolator reading interrupted!");
				Logger.errorException(ie);
			}
			System.out.println();
		}
	}

	public static double getFishersExactPvalue(float total1, float total2, float value1, float value2) {
		HypergeometricDistribution hyperDist=new HypergeometricDistribution(Math.round(total1+total2), Math.round(total1), Math.round(value1+value2));
		double pValue1=hyperDist.cumulativeProbability(Math.round(value1));

		hyperDist=new HypergeometricDistribution(Math.round(total1+total2), Math.round(total2), Math.round(value1+value2));
		double pValue2=hyperDist.cumulativeProbability(Math.round(value2));
		return Math.min(pValue1, pValue2);
	}

}
