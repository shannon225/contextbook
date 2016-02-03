package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import gnu.trove.list.array.TFloatArrayList;

public class PecanFeatureReaderTest {
	public static void main(String[] args) {
		File featureFile=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_980_1140.mzML.pecan.txt.features.txt");
		featureFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.mzML.pecan.txt.features.txt");
		//featureFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/yeast/Q_2014_0523_12_0_amol_uL_20mz.mzML.pecan.txt.features.txt");
		//featureFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/momo/20150301_A1_DIA_1.mzML.pecan.txt.features.txt");
		
		final HashMap<String, TFloatArrayList> targetData=new HashMap<String, TFloatArrayList>();
		final HashMap<String, TFloatArrayList> decoyData=new HashMap<String, TFloatArrayList>();
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				row.remove("id");
				row.remove("ScanNr");
				row.remove("topN");
				row.remove("precursorMz");
				row.remove("RTinMin");
				row.remove("sequence");
				row.remove("annotation");
				row.remove("charge2");
				row.remove("charge3");
				boolean isTarget=Integer.parseInt(row.remove("TD"))>0;
				
				for (Entry<String, String> entry : row.entrySet()) {
					String key=entry.getKey();
					try {
						float value=Float.parseFloat(entry.getValue());

						TFloatArrayList targets=targetData.get(key);
						TFloatArrayList decoys=decoyData.get(key);
						if (targets==null) {
							System.out.println("Got column: "+key);
							targets=new TFloatArrayList();
							targetData.put(key, targets);
							decoys=new TFloatArrayList();
							decoyData.put(key, decoys);
						}
						if (isTarget) {
							targets.add(value);
						} else {
							decoys.add(value);
						}
					} catch (NumberFormatException nfe) {
						System.err.println("error parsing ["+entry.getValue()+"] as a number from the ["+key+"] column!");
					}
				}
			}
		};

		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, featureFile, "\t", 1);
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

		TreeMap<String, ChartPanel> panelMap=new TreeMap<String, ChartPanel>();
		for (String key : targetData.keySet()) {
			TFloatArrayList targets=targetData.get(key);
			TFloatArrayList decoys=decoyData.get(key);
			
			ArrayList<XYPoint>[] points=PivotTableGenerator.createPivotTables(new float[][] {targets.toArray(), decoys.toArray()});
			XYTrace[] traces=new XYTrace[2];
			traces[0]=new XYTrace(points[0], GraphType.line, "Target");
			traces[1]=new XYTrace(points[1], GraphType.line, "Decoy");
			
			panelMap.put(key, Charter.getChart(key, "Count", true, traces));
		}
		Charter.launchCharts(featureFile.getName()+" Statistics", panelMap);
	}

}
