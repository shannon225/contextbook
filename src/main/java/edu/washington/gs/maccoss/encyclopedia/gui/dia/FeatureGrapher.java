package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import javax.swing.JFrame;
import javax.swing.JTabbedPane;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.Plot;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.FileChooserPanel;
import edu.washington.gs.maccoss.encyclopedia.gui.general.SimpleFilenameFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.PivotTableGenerator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest.RocPlot;
import edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest.RocPlot;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class FeatureGrapher {
	public static void main(String[] args) {
		File featureFile=FileChooserPanel.getFiles(null, "Feature text files", new SimpleFilenameFilter("features.txt"), (JFrame)null, true)[0];

		final JFrame f=new JFrame(featureFile.getName()+" Statistics");
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		f.getContentPane().add(graphFeatures(featureFile), BorderLayout.CENTER);

		f.pack();
		f.setSize(new Dimension(792, 612));
		f.setVisible(true);
	}

	public static JTabbedPane graphFeatures(File featureFile) {
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
				row.remove("Proteins");
				row.remove("charge1");
				row.remove("charge2");
				row.remove("charge3");
				row.remove("charge4");
				row.remove("protein");
				row.remove("pepLength");
				row.remove("averageFragmentDeltaMasses");
				row.remove("averageParentDeltaMass");
				row.remove("numConsidered");
				String targetDecoy = row.remove("TD");
				if (targetDecoy==null) {
					targetDecoy = row.remove("Label");
				}
				boolean isTarget=Integer.parseInt(targetDecoy)>0;
				
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

			@Override
			public void cleanup() {
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
		ArrayList<XYTrace> rocTraces=new ArrayList<>();
		for (String key : targetData.keySet()) {
			TFloatArrayList targets=targetData.get(key);
			TFloatArrayList decoys=decoyData.get(key);
			
			ArrayList<XYPoint>[] points=PivotTableGenerator.createPivotTables(new float[][] {targets.toArray(), decoys.toArray()}, true);
			XYTraceInterface[] traces=new XYTraceInterface[2];
			traces[0]=new XYTrace(points[0], GraphType.line, "Target");
			traces[1]=new XYTrace(points[1], GraphType.line, "Decoy");
			
			panelMap.put(key, Charter.getChart(key, "Count", true, traces));

			if (!key.equals("averageFragmentDeltaMasses")||!key.equals("averageParentDeltaMass")) {
				rocTraces.add(getRocPlot(key, targets.toArray(), decoys.toArray()));
			}
		}
		panelMap.put(" ROC", Charter.getChart("Q-Value", "Count", true, rocTraces.toArray(new XYTrace[rocTraces.size()])));
		return Charter.getTabbedChartPane(panelMap);
	}

	public static XYTrace getRocPlot(String name, float[] targets, float[] decoys) {
		if (name.equals("sumOfSquaredErrors")||name.equals("weightedSumOfSquaredErrors")||name.equals("precursorMass")||name.equals("percentBlankOverMono")||name.equals("numMissedCleavage")||name.equals("lnExpect")||name.equals("averageAbsParentDeltaMass")||name.equals("averageAbsFragmentDeltaMass")||name.equals("ms1MassError")||name.equals("ms2MassError")) {
			targets=General.multiply(targets, -1f);
			decoys=General.multiply(decoys, -1f);
		}

		ArrayList<ScoredObject<Boolean>> scores=new ArrayList<>();
		for (int i = 0; i < targets.length; i++) {
			scores.add(new ScoredObject<Boolean>(targets[i], true));
		}
		for (int i = 0; i < decoys.length; i++) {
			scores.add(new ScoredObject<Boolean>(decoys[i], false));
		}
		
		Collections.sort(scores);
		Collections.reverse(scores);

		int numTargets=0;
		int numDecoys=0;
		RocPlot plot=new RocPlot(name);
		for (ScoredObject<Boolean> score : scores) {
			if (score.y) {
				numTargets++;
			} else {
				numDecoys++;
			}

			float fdrValue=numTargets>0?numDecoys/(float)numTargets:1.0f;
			if (fdrValue>1.0f) fdrValue=1.0f;
			plot.addData(fdrValue, numTargets);
		}
		
		return plot.getTrace();
	}
}
