package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
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

public class StripeFileTrimmer {

	public static void main2(String[] args) throws Exception {
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);

		Range mzRange=new Range(399, 400);
		Range rtInSecRange=new Range(0, 2500000);
		File originalFile=new File("/Volumes/bcsbluessd/kkolotyuk/kkolotyuk_macos_run", "259070_GPF1.dia");
		File newFile=new File("/Volumes/bcsbluessd/kkolotyuk/kkolotyuk_macos_run", "trimmed_259070_GPF1.dia");
		trim(originalFile, newFile, mzRange, rtInSecRange, parameters);
	}

	public static void main(String[] args) throws Exception {
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);

		Range rtInSecRange=new Range(0, 2500000);
		File originalFile=new File("/Users/searleb/Documents/whoi/maggi", "191231_DIA_500to1100_metzyme_st5_80m.dia");

		chunk(originalFile, 10, rtInSecRange, parameters);
	}
	
	public static StripeFile trim(File originalFile, File newFile, Range mzRange, Range rtInSecRange, SearchParameters parameters) throws IOException, SQLException, DataFormatException {
		StripeFile stripeFile=new StripeFile(false);
		stripeFile.openFile();
		
		Logger.logLine("Splitting "+originalFile+" to selected from (mz:"+mzRange+" and rt:"+rtInSecRange+")...");
		StripeFileInterface thisStripeFile=StripeFileGenerator.getFile(originalFile, parameters);
		
		HashMap<Range, WindowData> dutyCycleMap=new HashMap<>();
		for (Entry<Range, WindowData> entry : thisStripeFile.getRanges().entrySet()) {
			if (mzRange.contains(entry.getKey().getMiddle())) {
				dutyCycleMap.put(entry.getKey(), entry.getValue());
			}
		}
		
		ArrayList<PrecursorScan> precursors = new ArrayList<>();
		for (PrecursorScan scan : thisStripeFile.getPrecursors(rtInSecRange.getStart(), rtInSecRange.getStop())) {
			precursors.add(scan);
		}
		stripeFile.addPrecursor(precursors);
		ArrayList<FragmentScan> stripes = new ArrayList<>();
		for (FragmentScan scan :  thisStripeFile.getStripes(mzRange, rtInSecRange.getStart(), rtInSecRange.getStop(), false)) {
			if (mzRange.contains(scan.getIsolationWindowLower())&&mzRange.contains(scan.getIsolationWindowUpper())) {
				stripes.add(scan);	
			}
		}
		stripeFile.addStripe(stripes);
		
		thisStripeFile.close();
		Logger.logLine("Finished merging, finalizing "+newFile.getName());

		stripeFile.setFileName(newFile.getName(), null, newFile.getAbsolutePath());
		stripeFile.setRanges(dutyCycleMap);

		stripeFile.saveAsFile(newFile);
		stripeFile.close();
		
		stripeFile=new StripeFile();
		stripeFile.openFile(newFile);
		return stripeFile;
	}
	public static void chunk(File originalFile, int numStripes, Range rtInSecRange, SearchParameters parameters) throws IOException, SQLException, DataFormatException {
		Logger.logLine("Splitting "+originalFile+" to break up into chunks of "+numStripes+" windows...");
		StripeFileInterface thisStripeFile=StripeFileGenerator.getFile(originalFile, parameters);
		
		Map<Range, WindowData> rangeMap = thisStripeFile.getRanges();
		ArrayList<Range> ranges=new ArrayList<>(rangeMap.keySet());
		Collections.sort(ranges);
		
		for (int i = 0; i < ranges.size(); i=i+numStripes) {		
			StripeFile stripeFile=new StripeFile(false);
			stripeFile.openFile();
			
			HashMap<Range, WindowData> dutyCycleMap=new HashMap<>();
			ArrayList<Range> targetRanges=new ArrayList<>();
			
			// find this chunks range
			double minMz=Double.MAX_VALUE;
			double maxMz=-Double.MAX_VALUE;
			for (int j = 0; j < numStripes; j++) {
				int index=i+j;
				if (i+j<ranges.size()) {
					Range range = ranges.get(index);
					targetRanges.add(range);
					dutyCycleMap.put(range, rangeMap.get(range));
					if (range.getStart()<minMz) {
						minMz=range.getStart();
					}
					if (range.getStop()>maxMz) {
						maxMz=range.getStop();
					}
				}
			}

			String absolutePath = originalFile.getAbsolutePath();
			absolutePath=absolutePath.substring(0, absolutePath.lastIndexOf('.'));
			File newFile=new File(absolutePath+"."+((int)Math.floor(minMz))+"."+((int)Math.floor(maxMz))+".dia");
			Logger.logLine("Working on "+newFile.getName());
			
			// add precursors that fit
			ArrayList<PrecursorScan> precursors = new ArrayList<>();
			for (PrecursorScan scan : thisStripeFile.getPrecursors(rtInSecRange.getStart(), rtInSecRange.getStop())) {
				Range precursorRange=new Range(scan.getIsolationWindowLower(), scan.getIsolationWindowUpper());
				if (precursorRange.contains(minMz)||precursorRange.contains(maxMz)) {
					precursors.add(scan);
				}
			}
			stripeFile.addPrecursor(precursors);
			
			for (Range mzRange : targetRanges) {
				ArrayList<FragmentScan> stripes = new ArrayList<>();
				for (FragmentScan scan :  thisStripeFile.getStripes(mzRange, rtInSecRange.getStart(), rtInSecRange.getStop(), false)) {
					if (mzRange.contains(scan.getIsolationWindowLower())&&mzRange.contains(scan.getIsolationWindowUpper())) {
						stripes.add(scan);	
					}
				}
				stripeFile.addStripe(stripes);
			}

			stripeFile.setFileName(newFile.getName(), null, newFile.getAbsolutePath());
			stripeFile.setRanges(dutyCycleMap);

			stripeFile.saveAsFile(newFile);
			stripeFile.close();
		}

		thisStripeFile.close();
	}
}
