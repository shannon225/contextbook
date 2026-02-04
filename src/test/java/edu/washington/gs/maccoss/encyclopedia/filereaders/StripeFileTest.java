package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.procedure.TIntFloatProcedure;

public class StripeFileTest {
	/** splits cycles into two files **/
	public static void main(String[] args) throws Exception {
		File diaFile=new File("/Users/searle.brian/Documents/temp/mapms_GPFDIA/raws/2024_06_17_125ng_HeLa_6xpGPFDIA_aurora_06.dia "); // here is the input raw file
		File f1File=new File("/Users/searle.brian/Documents/temp/mapms_GPFDIA/raws/2024_06_17_125ng_HeLa_6xpGPFDIA_aurora_06_first.dia "); 
		File f2File=new File("/Users/searle.brian/Documents/temp/mapms_GPFDIA/raws/2024_06_17_125ng_HeLa_6xpGPFDIA_aurora_06_second.dia "); 
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, PecanParameterParser.getDefaultParametersObject());
		
		ArrayList<PrecursorScan> precursors=stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE); // this reads all precursors with retention times from -INF to +INF
		Logger.logLine("Finished reading precursors");
		
		StripeFile firstCycle=new StripeFile(false);
		firstCycle.openFile();
		firstCycle.addPrecursor(precursors);
		HashMap<Range, TFloatArrayList> firstRetentionTimesByStripe=new HashMap<Range, TFloatArrayList>();

		StripeFile secondCycle=new StripeFile(false);
		secondCycle.openFile();
		secondCycle.addPrecursor(precursors);
		HashMap<Range, TFloatArrayList> secondRetentionTimesByStripe=new HashMap<Range, TFloatArrayList>();
		Logger.logLine("Finished writing precursors");
		
		ArrayList<FragmentScan> spectra=stripefile.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), -Float.MAX_VALUE, Float.MAX_VALUE, false);
		Logger.logLine("Finished reading MSMS");
		
		ArrayList<FragmentScan> thisCycle=new ArrayList<FragmentScan>();
		boolean firstFile=true;
		double lastTargetMz=0.0;
		for (FragmentScan msms : spectra) {
			double target=msms.getIsolationWindowCenter();
			if (target<lastTargetMz) { // cycle swap
				// publish cycle
				StripeFile file;
				HashMap<Range, TFloatArrayList> retentionTimesByStripe;
				if (firstFile) {
					file=firstCycle;
					retentionTimesByStripe=firstRetentionTimesByStripe;
					System.out.print("1");
				} else {
					file=secondCycle;
					retentionTimesByStripe=secondRetentionTimesByStripe;
					System.out.print("2");
				}
				file.addStripe(thisCycle);
				for (FragmentScan scan : thisCycle) {
					TFloatArrayList rtList=retentionTimesByStripe.get(scan.getRange());
					if (rtList==null) {
						rtList=new TFloatArrayList();
						retentionTimesByStripe.put(scan.getRange(), rtList);
					}
					rtList.add(scan.getScanStartTime());
				}
				
				firstFile=!firstFile;
				thisCycle.clear();
			}
			thisCycle.add(msms);
			lastTargetMz=target;
		}
		System.out.println();
		Logger.logLine("Finished writing MSMS");

		firstCycle.setFileName(f1File.getName(), null, f1File.getAbsolutePath());
		HashMap<Range, WindowData> firstDutyCycleMap=new HashMap<Range, WindowData>();
		for (Entry<Range, TFloatArrayList> entry : firstRetentionTimesByStripe.entrySet()) {
			Range range=entry.getKey();
			TFloatArrayList rts=entry.getValue();
			float[] deltas=General.firstDerivative(rts.toArray());
			float averageDutyCycle=General.mean(deltas);
			firstDutyCycleMap.put(range, new WindowData(averageDutyCycle, rts.size()));
		}
		firstCycle.setRanges(firstDutyCycleMap);

		firstCycle.saveAsFile(f1File);
		firstCycle.close();

		secondCycle.setFileName(f2File.getName(), null, f2File.getAbsolutePath());
		HashMap<Range, WindowData> secondDutyCycleMap=new HashMap<Range, WindowData>();
		for (Entry<Range, TFloatArrayList> entry : secondRetentionTimesByStripe.entrySet()) {
			Range range=entry.getKey();
			TFloatArrayList rts=entry.getValue();
			float[] deltas=General.firstDerivative(rts.toArray());
			float averageDutyCycle=General.mean(deltas);
			secondDutyCycleMap.put(range, new WindowData(averageDutyCycle, rts.size()));
		}
		secondCycle.setRanges(secondDutyCycleMap);

		secondCycle.saveAsFile(f2File);
		secondCycle.close();

		Logger.logLine("Done!");
	}
	
	
	public static void main2(String[] args) throws Exception {
		File diaFile=new File("/Users/searle.brian/Downloads/2024_06_17_125ng_HeLa_1e6_6xpDIA_aurora_01.dia"); // here is the input raw file
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, PecanParameterParser.getDefaultParametersObject());
		
		ArrayList<PrecursorScan> scans=stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE); // this reads all precursors with retention times from -INF to +INF
		ArrayList<XYPoint> tic=new ArrayList<XYPoint>();
		TIntFloatHashMap mzSum=new TIntFloatHashMap();
		for (PrecursorScan scan : scans) { // for each MS1 (PrecursorScan)
			float rt=scan.getScanStartTime()/60f; // get the retention time in minutes
			double[] masses=scan.getMassArray(); // the mass and intensity arrays have the same length (e.g., mass[i] corresponds to intensity[i])
			float[] intensities=scan.getIntensityArray();
			float intensity=General.sum(intensities);
			XYPoint point=new XYPoint(rt, intensity);
			tic.add(point);
			
			for (int i=0; i<intensities.length; i++) {
				int index=(int)Math.round(masses[i]);
				mzSum.adjustOrPutValue(index, intensities[i], intensities[i]);
			}
		}

		final ArrayList<XYPoint> mzs=new ArrayList<XYPoint>();
		mzSum.forEachEntry(new TIntFloatProcedure() {
			@Override
			public boolean execute(int arg0, float arg1) {
				mzs.add(new XYPoint(arg0, arg1));
				System.out.println(arg0+"\t"+arg1);
				return true;
			}
		});

		//Charter.launchChart("Retention Time (Minutes)", "Total Ion Chromatogram", false, new XYTrace(tic, GraphType.line, "TIC"));
		Charter.launchChart("M/z", "Total Ion Chromatogram", false, new XYTrace(mzs, GraphType.line, "Mz"));
	}
}

