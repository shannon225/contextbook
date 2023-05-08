package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.math.MedianInterpolatorTest;
import junit.framework.TestCase;

public class MassErrorKDETest extends TestCase {

	public void testErrors() {
		InputStream is=MedianInterpolatorTest.class.getResourceAsStream("/deltamz/mass_errors_example.txt");
		
		final ArrayList<XYPoint> fragments=new ArrayList<XYPoint>();
		final ArrayList<XYPoint> parents=new ArrayList<XYPoint>();

		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				float rt=Float.parseFloat(row.get("RTinMin"));
				fragments.add(new XYPoint(rt, Float.parseFloat(row.get("averageFragmentDeltaMasses"))));
				parents.add(new XYPoint(rt, Float.parseFloat(row.get("averageParentDeltaMass"))));
			}
			
			@Override
			public void cleanup() {
			}
		};

		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, is, "\t", 1);
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

		long time=System.currentTimeMillis();
		MassErrorFilter filter=MassErrorFilter.getFilter(new MassTolerance(10), 1, fragments);
		System.out.println((System.currentTimeMillis()-time)+" ms");
		
		Range rt=new Range(XYTrace.toFloatArrays(parents).x);
		assertEquals(1, rt.getStart(), 1);
		assertEquals(115, rt.getStop(), 1);
		
		for (int i = Math.round(rt.getStart()); i <= Math.round(rt.getStop()); i++) {
			assertEquals(0.943529521f, filter.getYValue(i), 0.5f);
		}

		time=System.currentTimeMillis();
		filter=MassErrorFilter.getFilter(new MassTolerance(10), 1, parents);
		System.out.println((System.currentTimeMillis()-time)+" ms");
		
		for (int i = Math.round(rt.getStart()); i <= Math.round(rt.getStop()); i++) {
			assertEquals(2.892641459, filter.getYValue(i), 1f);
		}
	}

}
