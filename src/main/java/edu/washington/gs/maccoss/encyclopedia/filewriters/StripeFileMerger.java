package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Set;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;

public class StripeFileMerger {
	private static final String FILENAME_DELIMITER = ";";
	
	public static void main(String[] args) throws Exception {
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);
		
//		File dir=new File("/Users/searleb/Downloads/wide_iso");
//		
//		merge(new File[] {new File(dir, "Loo_2020_0406_RJ_102_32.mzML"), new File(dir,"Loo_2020_0406_RJ_103_32.mzML")}, 
//				new File(dir, "combined.dia"), parameters);
		

		File dir=new File("/Volumes/searle_ssd/cobbs/SN606_RIPA/");
		
		merge(new File[] {new File(dir, "20200708_BCS_LOOM_COBBS_GPFDIA_SN606_1.mzML"), new File(dir,"20200708_BCS_LOOM_COBBS_GPFDIA_SN606_2.mzML"),
				new File(dir,"20200708_BCS_LOOM_COBBS_GPFDIA_SN606_3.mzML"), new File(dir,"20200708_BCS_LOOM_COBBS_GPFDIA_SN606_4.mzML"),
				new File(dir,"20200708_BCS_LOOM_COBBS_GPFDIA_SN606_5.mzML"), new File(dir,"20200708_BCS_LOOM_COBBS_GPFDIA_SN606_6.mzML")}, 
				new File(dir, "SN606_combined.dia"), parameters);
	}
	public static StripeFile merge(File[] fs, File newFile, SearchParameters parameters) throws IOException, SQLException, DataFormatException {
		StripeFile stripeFile=new StripeFile();
		stripeFile.openFile();
		HashMap<Range, WindowData> dutyCycleMap=new HashMap<>();
		
		int scanIndex=0;
		for (int i = 0; i < fs.length; i++) {
			Logger.logLine("Adding "+fs[i].getName()+" to merged file ("+(i+1)+" of "+fs.length+")...");
			StripeFileInterface thisStripeFile=StripeFileGenerator.getFile(fs[i], parameters);
			Map<Range, WindowData> ranges = thisStripeFile.getRanges();
			TFloatArrayList timeBetweenScans=new TFloatArrayList();
			for (WindowData time : ranges.values()) {
				timeBetweenScans.add(time.getAverageDutyCycle());
			}
			float maxTimeForCommon=1.5f*QuickMedian.median(timeBetweenScans.toArray());
			HashMap<Range, WindowData> primaryRanges = new HashMap<Range, WindowData>(ranges);
			for (Entry<Range, WindowData> entry : ranges.entrySet()) {
				if (entry.getValue().getAverageDutyCycle()>maxTimeForCommon) {
					primaryRanges.remove(entry.getKey());
				}
			}
			Range mzRange=Range.getWidestRange(new ArrayList<Range>(primaryRanges.keySet()));
			Logger.logLine("Keeping m/z range of "+mzRange.toString()+"...");
			
			dutyCycleMap.putAll(primaryRanges);
			ArrayList<PrecursorScan> precursors = new ArrayList<>();
			for (PrecursorScan scan : thisStripeFile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE)) {
				precursors.add(scan.shallowClone(i, scanIndex++, mzRange));
			}
			stripeFile.addPrecursor(precursors);
			ArrayList<FragmentScan> stripes = new ArrayList<>();
			for (FragmentScan scan :  thisStripeFile.getStripes(mzRange, -Float.MAX_VALUE, Float.MAX_VALUE, false)) {
				stripes.add(scan.shallowClone(i, scanIndex++));
			}
			stripeFile.addStripe(stripes);
			
			thisStripeFile.close();
		}
		Logger.logLine("Finished merging, finalizing "+newFile.getName());

		stripeFile.setFileName(newFile.getName(), null, newFile.getAbsolutePath());
		stripeFile.setRanges(dutyCycleMap);

		stripeFile.saveAsFile(newFile);
		stripeFile.close();
		
		stripeFile=new StripeFile();
		stripeFile.openFile(newFile);
		return stripeFile;
	}
}
