package edu.washington.gs.maccoss.encyclopedia.gui.dia;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.StringTokenizer;
import java.util.TreeMap;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;

import edu.washington.gs.maccoss.encyclopedia.datastructures.GlobalRangeTracker;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ScanRangeTracker;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlScanRangeTrackerSAXProducer;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.gui.general.GUIParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;

public class MzmlStructureCharter {

	private static final int MAXIMUM_NUMBER_OF_SCANS_PER_TYPE = 500;
	
	public static void main(String[] args) throws Exception {
		File f=new File("/Users/searleb/Downloads/Mass List Table3.csv");
		ScanRangeTracker tracker=new ScanRangeTracker();
		

		BufferedReader in=new BufferedReader(new FileReader(f));
		String eachline;
		int count=-1;
		while ((eachline=in.readLine())!=null) {
			count++;
			if (count<=0||eachline.trim().length()==0) {
				continue;
			}
			
			StringTokenizer st=new StringTokenizer(eachline, "-");
			double start=Double.parseDouble(st.nextToken());
			double stop=Double.parseDouble(st.nextToken());
			System.out.println(count+") "+start+" to "+stop);
			tracker.addRange(new Range(start, stop), count);
		}
		in.close();
		ChartPanel panel=MzmlStructureCharter.getStructureChart(tracker, true);
		Charter.launchComponent(panel, "File structure", new Dimension(900, 450));
	}

	public static ChartPanel getGlobalStructureChart(StripeFile dia) {
		try {
			Connection c=dia.getConnection();
			try {
				Statement s=c.createStatement();
				try {
					GlobalRangeTracker tracker=new GlobalRangeTracker();
					Logger.logLine("Strarting to read windows...");
					
					// double[] is {stopMz, minRT, maxRT}
					TDoubleObjectHashMap<double[]> valuesByLowerBound=new TDoubleObjectHashMap<double[]>();

					ResultSet rs=s.executeQuery("select scanstarttime, isolationwindowlower, isolationwindowupper from spectra");
					while (rs.next()) {
						double rt=rs.getFloat(1);
						double startMz=rs.getDouble(2);
						double stopMz=rs.getDouble(3);
						
						double[] values=valuesByLowerBound.get(startMz);
						if (values==null) {
							values=new double[] {stopMz, rt, rt};
							valuesByLowerBound.put(startMz, values);
						} else {
							if (values[1]>rt) values[1]=rt;
							if (values[2]<rt) values[2]=rt;
						}
					}
					
					valuesByLowerBound.forEachEntry(new TDoubleObjectProcedure<double[]>() {
						@Override
						public boolean execute(double a, double[] b) {
							tracker.addRange(new Range(a, b[0]), new Range(b[1], b[2]));
							return true;
						}
					});
					
					Logger.logLine("Found "+tracker.getStripeRTsInSecs().size()+" total windows...");
					
					return getStructureChart(tracker, false);

				} finally {
					s.close();
				}
			} finally {
				c.close();
			}
			
		} catch (IOException ioe) {
			throw new EncyclopediaException("DIA reading IO error!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("DIA reading SQL error!", sqle);
		}
	}
	
	public static ChartPanel getStructureChart(StripeFile dia) {
		try {
			Connection c=dia.getConnection();
			try {
				Statement s=c.createStatement();
				try {
					float maxRT=Float.MAX_VALUE;
					ScanRangeTracker tracker=new ScanRangeTracker();
					double maxUpper=0.0;
					double minLower=Float.MAX_VALUE;
					
					ResultSet rs=s.executeQuery("select scanstarttime, isolationwindowlower, isolationwindowupper from spectra order by scanstarttime limit 1000");
					while (rs.next()) {
						float scanStartTime=rs.getFloat(1);
						double isolationWindowLower=rs.getDouble(2);
						double isolationWindowUpper=rs.getDouble(3);
						if (isolationWindowLower<minLower) minLower=isolationWindowLower;
						if (isolationWindowUpper>maxUpper) maxUpper=isolationWindowUpper;
						
						boolean keepGoing=tracker.addRange(new Range(isolationWindowLower, isolationWindowUpper), scanStartTime);
						if (!keepGoing) {
							maxRT=scanStartTime;
							break;
						}
					}
					

					rs=s.executeQuery("select scanstarttime, isolationwindowlower, isolationwindowupper from precursor where scanstarttime<="+maxRT+" order by scanstarttime");
					while (rs.next()) {
						float scanStartTime=rs.getFloat(1);
						double isolationWindowLower=rs.getDouble(2);
						double isolationWindowUpper=rs.getDouble(3);
						if (isolationWindowUpper>1e8) { 
							isolationWindowLower=minLower;
							isolationWindowUpper=maxUpper;
						}
						tracker.addPrecursor(new Range(isolationWindowLower, isolationWindowUpper), scanStartTime);
					}
					
					return getStructureChart(tracker, false);

				} finally {
					s.close();
				}
			} finally {
				c.close();
			}
			
		} catch (IOException ioe) {
			throw new EncyclopediaException("DIA reading IO error!", ioe);
		} catch (SQLException sqle) {
			throw new EncyclopediaException("DIA reading SQL error!", sqle);
		}
	}

	public static ChartPanel getStructureChart(File mzMLFile, boolean isGlobal) {
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);

		if (mzMLFile.getName().toLowerCase().endsWith("dia")) {
			StripeFileInterface dia=StripeFileGenerator.getFile(mzMLFile, parameters);
			if (dia instanceof StripeFile) {

				if (isGlobal) {
					return getGlobalStructureChart((StripeFile)dia);
				} else {
					return getStructureChart((StripeFile)dia);
				}
			}
		}

		if (isGlobal) {
			Logger.errorLine("Missing structure chart because file is not already built.");
			ChartPanel panel=Charter.getChart("M/Z", "Retention Time", false, new XYTraceInterface[0]);
			return panel;
		}

		ScanRangeTracker scanTracker=null;
		Logger.logLine("Indexing "+mzMLFile.getName()+" ...");

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

		return getStructureChart(scanTracker, false);
	}

	public static ChartPanel getStructureChart(ScanRangeTracker scanTracker, boolean isScanNumberInsteadOfRT) {
		TreeMap<Range, TFloatArrayList> retentionTimesByStripe=new TreeMap<>(scanTracker.getStripeRTsInSecs());
		TreeMap<Range, TFloatArrayList> retentionTimesByPrecursor=new TreeMap<>(scanTracker.getPrecursorRTsInSecs());

		float firstScan=Float.MAX_VALUE;
		float lastScan=0.0f;
		ArrayList<XYTraceInterface> traces=new ArrayList<>();
		boolean everyOther=false;
		int totalAllowableFragments=MAXIMUM_NUMBER_OF_SCANS_PER_TYPE;
		for (Entry<Range, TFloatArrayList> entry : retentionTimesByStripe.entrySet()) {
			Range range=entry.getKey();
			TFloatArrayList rts=entry.getValue();
			if (rts.size()>0) {
				totalAllowableFragments--;
				if (totalAllowableFragments<=0) break;
				
				float rt = rts.get(0);
				if (rt<firstScan) firstScan=rt;
				if (rt>lastScan) lastScan=rt;
				everyOther=!everyOther;

				XYTraceInterface trace=new XYTrace(new float[] {range.getStart(), range.getStop()}, new float[] {rt, rt}, GraphType.squaredline, range.toString(), getColor(everyOther), 5.0f);
				traces.add(trace);
				if (rts.size()>1) {
					float secondRT = rts.get(1);
					if (secondRT<firstScan) firstScan=secondRT;
					if (secondRT>lastScan) lastScan=secondRT;
					trace=new XYTrace(new float[] {range.getStart(), range.getStop()}, new float[] {secondRT, secondRT}, GraphType.squaredline, range.toString(), getColor(everyOther), 5.0f);
					traces.add(trace);
					trace=new XYTrace(new float[] {range.getStop(), range.getStop()}, new float[] {rt, secondRT}, GraphType.dashedline, range.toString(), Color.gray, 1.0f);
					traces.add(trace);
				}
			}
		}
		
		float rtRangeMargin=(lastScan-firstScan)*0.2f;
		Range rtRange=new Range(firstScan-rtRangeMargin, lastScan+rtRangeMargin);

		int totalAllowablePrecursors=MAXIMUM_NUMBER_OF_SCANS_PER_TYPE;
		for (Entry<Range, TFloatArrayList> entry : retentionTimesByPrecursor.entrySet()) {
			Range range=entry.getKey();
			TFloatArrayList rts=entry.getValue();
			for (float rt : rts.toArray()) {
				if (rtRange.contains(rt)) {
					totalAllowablePrecursors--;
					if (totalAllowablePrecursors<=0) break;
					
					XYTraceInterface trace=new XYTrace(new float[] {range.getStart(), range.getStop()}, new float[] {rt, rt}, GraphType.squaredline, range.toString(), Color.LIGHT_GRAY, 5.0f);
					traces.add(trace);
				}
			}
		}
		Logger.errorLine("FOUND "+traces.size()+" TRACES");

		String yAxis=isScanNumberInsteadOfRT?"Scan Number":"Retention Time (secs)";
		ChartPanel panel=Charter.getChart("M/Z", yAxis, false, traces.toArray(new XYTraceInterface[traces.size()]));
		return panel;
	}


	public static ChartPanel getStructureChart(GlobalRangeTracker scanTracker, boolean isScanNumberInsteadOfRT) {
		TreeMap<Range, Range> retentionTimesByStripe=new TreeMap<>(scanTracker.getStripeRTsInSecs());
		if (retentionTimesByStripe.size()==0) {
			// FIXME should include precursor data eventually
			return null;
		}

		boolean everyOther=false;

		String yAxis=isScanNumberInsteadOfRT?"Scan Number":"Retention Time (min)";

		ArrayList<XYShapeAnnotation> shapes=new ArrayList<XYShapeAnnotation>();
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		double minMz=Double.MAX_VALUE;
		double minRT=Double.MAX_VALUE;
		double maxMz=Double.MIN_VALUE;
		double maxRT=Double.MIN_VALUE;
		
		for (Entry<Range, Range> entry : retentionTimesByStripe.entrySet()) {
			Range mzRange=entry.getKey();
			Range rtRange=entry.getValue();
			if (mzRange.getStart()<minMz) minMz=mzRange.getStart();
			if (mzRange.getStop()>maxMz) maxMz=mzRange.getStop();
			if (rtRange.getStart()<minRT) minRT=rtRange.getStart();
			if (rtRange.getStop()>maxRT) maxRT=rtRange.getStop();
			
			double x=mzRange.getStart();
			double y=rtRange.getStart()/60f;
			double width=mzRange.getStop()-mzRange.getStart();
			double height=(rtRange.getStop()-rtRange.getStart())/60f;
			Rectangle2D shape = new Rectangle2D.Double();
			shape.setFrame(x, y, width, height);
			
			everyOther=!everyOther;
			
			BasicStroke stroke=new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
			shapes.add(new XYShapeAnnotation(shape, stroke, Color.gray, getColor50p(everyOther)));
		}
		
		points.add(new XYPoint(minMz, minRT/60.0));
		points.add(new XYPoint(maxMz, minRT/60.0));
		points.add(new XYPoint(maxMz, maxRT/60.0));
		points.add(new XYPoint(minMz, maxRT/60.0));
		points.add(new XYPoint(minMz, minRT/60.0));
		
		return Charter.getShapeChart(null, "M/Z", yAxis, 16, 16, shapes, points, false);
	}

	private static Color getColor(boolean everyOther) {
		//return everyOther?new Color(0, 0, 200):new Color(100, 100, 255);
		return everyOther?GUIParameters.getBaseColor():GUIParameters.getBrighterColor();
	}

	private static Color getColor50p(boolean everyOther) {
		//return everyOther?new Color(0, 0, 200):new Color(100, 100, 255);
		return everyOther?GUIParameters.getBaseColor(127):GUIParameters.getBrighterColor(127);
	}
}
