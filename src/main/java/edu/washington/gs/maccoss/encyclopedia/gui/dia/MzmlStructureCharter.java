package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ScanRangeTracker;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlScanRangeTrackerSAXProducer;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import gnu.trove.list.array.TFloatArrayList;

public class MzmlStructureCharter {

	public static ChartPanel getStructureChart(File mzMLFile) {
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);

		ScanRangeTracker scanTracker=null;
		try {
			Logger.logLine("Indexing "+mzMLFile.getName()+" ...");
			StripeFile stripeFile=new StripeFile();
			stripeFile.openFile();

			MzmlScanRangeTrackerSAXProducer producer=new MzmlScanRangeTrackerSAXProducer(mzMLFile, parameters);
			scanTracker=producer.getRetentionTimesByStripe();

			Thread producerThread=new Thread(producer);

			Thread[] threads=new Thread[] {producerThread};

			for (int i=0; i<threads.length; i++) {
				threads[i].start();
			}

			try {
				for (int i=0; i<threads.length; i++) {
					threads[i].join();
				}
			} catch (InterruptedException ie) {
				Logger.errorLine("DIA reading interrupted!");
				Logger.errorException(ie);
			}
		} catch (IOException ioe) {
			throw new EncyclopediaException("DIA reading IO error!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("DIA reading SQL error!", sqle);
		}
		
		return getStructureChart(scanTracker, false);
	}

	public static ChartPanel getStructureChart(ScanRangeTracker scanTracker, boolean isScanNumberInsteadOfRT) {
		TreeMap<Range, TFloatArrayList> retentionTimesByStripe=new TreeMap<>(scanTracker.getStripeRTsInSecs());
		TreeMap<Range, TFloatArrayList> retentionTimesByPrecursor=new TreeMap<>(scanTracker.getPrecursorRTsInSecs());

		float firstScan=Float.MAX_VALUE;
		ArrayList<XYTraceInterface> traces=new ArrayList<>();
		boolean everyOther=false;
		for (Entry<Range, TFloatArrayList> entry : retentionTimesByStripe.entrySet()) {
			Range range=entry.getKey();
			TFloatArrayList rts=entry.getValue();
			if (rts.size()>0) {
				if (rts.get(0)<firstScan) firstScan=rts.get(0);
				everyOther=!everyOther;

				XYTraceInterface trace=new XYTrace(new float[] {range.getStart(), range.getStop()}, new float[] {rts.get(0), rts.get(0)}, GraphType.squaredline, range.toString(), getColor(everyOther), 5.0f);
				traces.add(trace);
				if (rts.size()>1) {
					trace=new XYTrace(new float[] {range.getStart(), range.getStop()}, new float[] {rts.get(1), rts.get(1)}, GraphType.squaredline, range.toString(), getColor(everyOther), 5.0f);
					traces.add(trace);
					trace=new XYTrace(new float[] {range.getStop(), range.getStop()}, new float[] {rts.get(0), rts.get(1)}, GraphType.dashedline, range.toString(), Color.gray, 1.0f);
					traces.add(trace);
				}
			}
		}

		for (Entry<Range, TFloatArrayList> entry : retentionTimesByPrecursor.entrySet()) {
			Range range=entry.getKey();
			TFloatArrayList rts=entry.getValue();
			for (float rt : rts.toArray()) {
				if (rt>firstScan) {
					XYTraceInterface trace=new XYTrace(new float[] {range.getStart(), range.getStop()}, new float[] {rt, rt}, GraphType.squaredline, range.toString(), Color.LIGHT_GRAY, 5.0f);
					traces.add(trace);
				}
			}
		}

		String yAxis=isScanNumberInsteadOfRT?"Scan Number":"Retention Time (secs)";
		ChartPanel panel=Charter.getChart("M/Z", yAxis, false, traces.toArray(new XYTraceInterface[traces.size()]));
		return panel;
	}

	private static Color getColor(boolean everyOther) {
		return everyOther?new Color(0, 0, 200):new Color(100, 100, 255);
	}
}
