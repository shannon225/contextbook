package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
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

public class KGGAlignmentTest {
	public static void main(String[] args) throws Exception {
		ArrayList<XYPoint> rts=getKGGData();
		//ArrayList<XYPoint> rts=getNonKGGData();
		RetentionTimeFilter filter=RetentionTimeFilter.getFilter(rts);
		filter.plot(rts, Optional.ofNullable((File)null));
	}

	public static ArrayList<XYPoint> getKGGData() throws Exception {
		File file=new File("/Users/searleb/Documents/sperm/ubiquitin/U2SO-Untreated_RT_KGG_DBW210414-2226.tsv");
		return getData(new FileInputStream(file));
	}
	public static ArrayList<XYPoint> getNonKGGData() throws Exception {
		File file=new File("/Users/searleb/Documents/sperm/ubiquitin/U2SO-Untreated_RT_nonKGG_DBW210414-2226.tsv");
		return getData(new FileInputStream(file));
	}

	public static ArrayList<XYPoint> getData(InputStream is) {
		final ArrayList<XYPoint> rts=new ArrayList<XYPoint>();

		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				float predicted=Float.parseFloat(row.get("Prosit_RT"));
				float actual=Float.parseFloat(row.get("DDA_RT"));
				rts.add(new XYPoint(predicted/60f, actual/60f));
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

		return rts;
	}
}
