package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import gnu.trove.list.array.TFloatArrayList;

public class MzmlToDIAProducerTest {
	
	public static void main(String[] args) {
		File mzMLFile=new File("/Users/searleb/Documents/school/projects/pecandata/121115_BCS_HeLa_24mz_400_1000.mzML");

		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);

		try {
			Logger.logLine("Indexing "+mzMLFile.getName()+" ...");
			StripeFile stripeFile=new StripeFile();
			stripeFile.openFile();

			BlockingQueue<MzmlBlock> mzmlBlockQueue=new LinkedBlockingQueue<MzmlBlock>();
			MzmlToDIASAXProducer producer=new MzmlToDIASAXProducer(mzMLFile, mzmlBlockQueue, parameters);
			HashMap<Range, TFloatArrayList> retentionTimesByStripe=producer.getRetentionTimesByStripe();
			
			Thread producerThread=new Thread(producer);

			Thread[] threads=new Thread[] { producerThread };

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
			
			ArrayList<XYTraceInterface> traces=new ArrayList<>();
			for (Entry<Range, TFloatArrayList> entry : retentionTimesByStripe.entrySet()) {
				Range range=entry.getKey();
				TFloatArrayList rts=entry.getValue();
				if (rts.size()>0) {
					XYTraceInterface trace=new XYTrace(new float[] { range.getStart(), range.getStop() }, new float[] { rts.get(0), rts.get(0) }, GraphType.boldline, range.toString(), Color.BLUE, 5.0f);
					traces.add(trace);
					if (rts.size()>1) {
						trace=new XYTrace(new float[] { range.getStart(), range.getStop() }, new float[] { rts.get(1), rts.get(1) }, GraphType.boldline, range.toString(), Color.blue.darker(), 5.0f);	
						traces.add(trace);
						trace=new XYTrace(new float[] { range.getStop(), range.getStop() }, new float[] { rts.get(0), rts.get(1) }, GraphType.dashedline, range.toString(), Color.gray, 2.0f);
						traces.add(trace);
					}
				}
			}
			
			Charter.launchChart("M/Z", "Retention Time (secs)", false, traces.toArray(new XYTraceInterface[traces.size()]));

		} catch (IOException ioe) {
			throw new EncyclopediaException("DIA reading IO error!", ioe);
		} catch (SQLException sqle) {
			sqle.printStackTrace();
			throw new EncyclopediaException("DIA reading SQL error!", sqle);
		}
	}
}
