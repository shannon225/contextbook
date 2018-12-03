package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.InputStream;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.gui.dia.WindowingSchemeWizard;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Histogram;
import gnu.trove.list.array.TFloatArrayList;

public class WindowSchemeGenerator {
	public static final int NORMAL_DIA=0; 
	public static final int OVERLAP_DIA=1; 
	public static final int VARIABLE_WIDTH_DIA=2; 
	public static final int VARIABLE_WIDTH_PHOSPHO_DIA=3;
	
	private static final double optimalMzConstant=0.25;
	private static final double optimalMzIncrement=1.00045475;
	private static Histogram normalProteome=null;
	private static Histogram phosphoProteome=null; 
	
	/**
	 * 
	 * @param windowingSchemeIndex
	 * @param start
	 * @param stop
	 * @param numWindows
	 * @param margin
	 * @return
	 */
	public static ScanRangeTracker generateWindowingScheme(int windowingSchemeIndex, int start, int stop, int numWindows, float margin) {
		ScanRangeTracker tracker=new ScanRangeTracker();
		if (windowingSchemeIndex==NORMAL_DIA) {
			// normal DIA
			int increment=(int)Math.ceil((stop-start)/(float)numWindows);
			for (int i=0; i<numWindows; i++) {
				float left=getOptimalBoundary(start+increment*i)-margin;
				float right=getOptimalBoundary(start+increment*(i+1))+margin;
				tracker.addRange(new Range(left, right), i+1);
			}
		} else if (windowingSchemeIndex==OVERLAP_DIA) {
			// overlap DIA
			int increment=(int)Math.ceil((stop-start)/(float)numWindows);
			for (int i=0; i<numWindows; i++) {
				float left=getOptimalBoundary(start+increment*i)-margin;
				float right=getOptimalBoundary(start+increment*(i+1))+margin;
				tracker.addRange(new Range(left, right), i+1);
			}
			for (int i=0; i<numWindows-1; i++) {
				float left=getOptimalBoundary(start+increment*(i+0.5f))-margin;
				float right=getOptimalBoundary(start+increment*(i+1.5f))+margin;
				tracker.addRange(new Range(left, right), numWindows+i+1);
			}
		} else if (windowingSchemeIndex==VARIABLE_WIDTH_DIA) {
			// variable width normal DIA
			Histogram hist=getNormalProteome();
			trackFromHistogram(start, stop, numWindows, margin, tracker, hist);
			
		} else if (windowingSchemeIndex==VARIABLE_WIDTH_PHOSPHO_DIA) {
			// variable width phospho DIA
			Histogram hist=getPhosphoProteome();
			trackFromHistogram(start, stop, numWindows, margin, tracker, hist);
		}
		return tracker;
	}

	private static void trackFromHistogram(int start, int stop, int numWindows, float margin, ScanRangeTracker tracker, Histogram hist) {
		float[] targetPercentiles=new float[numWindows+1];
		for (int i=0; i<targetPercentiles.length; i++) {
			targetPercentiles[i]=i/(float)numWindows;
		}
		targetPercentiles[numWindows]=1.0f;
		
		float[] boundaries=hist.getPercentiles(targetPercentiles, start, stop);
		for (int i=1; i<boundaries.length; i++) {
			float left=getOptimalBoundary(boundaries[i-1])-margin;
			float right=getOptimalBoundary(boundaries[i])+margin;
			tracker.addRange(new Range(left, right), i);
		}
	}
	
	private static float getOptimalBoundary(float nominalMass) {
		return (float)(Math.ceil(nominalMass/optimalMzIncrement)*optimalMzIncrement+optimalMzConstant);
	}

	public static Histogram getPhosphoProteome() {
		if (phosphoProteome==null) {
			InputStream is=WindowingSchemeWizard.class.getResourceAsStream("/data/hela_phospho_mzs.txt");
			phosphoProteome=getData(is);
		}
		return phosphoProteome;
	}
	public static Histogram getNormalProteome() {
		if (normalProteome==null) {
			InputStream is=WindowingSchemeWizard.class.getResourceAsStream("/data/hela_mzs.txt");
			normalProteome=getData(is);
		}
		return normalProteome;
	}

	public static Histogram getData(InputStream is) {
		TFloatArrayList mzs=new TFloatArrayList();
		TFloatArrayList counts=new TFloatArrayList();

		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				float mz=Float.parseFloat(row.get("mz"));
				float count=Float.parseFloat(row.get("count"));
				mzs.add(mz);
				counts.add(count);
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
			Logger.errorLine("Histogram reading interrupted!");
			Logger.errorException(ie);
		}

		return Histogram.generateHistogram(mzs.toArray(), counts.toArray(), 1.0f);
	}
}
