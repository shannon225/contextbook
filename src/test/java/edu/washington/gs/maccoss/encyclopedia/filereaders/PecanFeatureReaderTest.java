package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import javax.swing.JFrame;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
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
		File featureFile=FileChooserPanel.getFiles(null, "Feature text files", new SimpleFilenameFilter("features.txt"), (JFrame)null)[0];
		
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
