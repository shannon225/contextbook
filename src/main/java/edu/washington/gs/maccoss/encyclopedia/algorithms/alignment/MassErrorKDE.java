package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;

public class MassErrorKDE extends RectangleKDE {
	public static final int MASS_ERROR_RESOLUTION = 1000;
	public static final int RT_RESOLUTION = 3000;

	public MassErrorKDE(ArrayList<XYPoint> points) {
		super(points, RT_RESOLUTION, MASS_ERROR_RESOLUTION);
	}
	
	public static void main(String[] args) {
		File f=new File("/Users/searleb/Documents/encyclopedia/mass_errors_example.txt");

		final ArrayList<XYPoint> targets=new ArrayList<XYPoint>();
		final ArrayList<XYPoint> decoys=new ArrayList<XYPoint>();

		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				int label=Integer.parseInt(row.get("Label"));
				float rt=Float.parseFloat(row.get("RTinMin"));
				//float delta=Float.parseFloat(row.get("averageFragmentDeltaMasses"));
				float delta=Float.parseFloat(row.get("averageParentDeltaMass"));
				if (delta!=10) {
					if (label==1) {
						targets.add(new XYPoint(rt, delta));
					} else {
						decoys.add(new XYPoint(rt, delta));
					}
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

		MassErrorFilter filter=MassErrorFilter.getFilter(new MassTolerance(10), 1, targets);
		filter.plot(targets, Optional.empty());

//		ArrayList<XYTrace> traces=new ArrayList<>();
//		MassErrorKDE filter=null;
//		for (int i=0; i<1; i++) {
//			long time=System.currentTimeMillis();
//			filter=new MassErrorKDE(targets);
//			Function func=filter.trace();
//			XYTrace knots=new XYTrace(func.getKnots(), GraphType.line, "Iter:"+i);
//			traces.add(knots);
//			System.out.println((System.currentTimeMillis()-time));
//		}
//		filter.plot();
//
//		ArrayList<XYTrace> targetTraces=new ArrayList<>(traces);
//		targetTraces.add(new XYTrace(targets, GraphType.tinypoint, "Targets", new Color(0, 0, 255, 20), 1.0f));
//		Charter.launchChart("RT", "PPM", false, targetTraces.toArray(new XYTrace[0]));
//
//		ArrayList<XYTrace> decoyTraces=new ArrayList<>(traces);
//		decoyTraces.add(new XYTrace(decoys, GraphType.tinypoint, "Decoys", new Color(255, 0, 0, 20), 1.0f));
//		Charter.launchChart("RT", "PPM", false, decoyTraces.toArray(new XYTrace[0]));
	}
}
