package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
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

	public static void main(String[] args) throws Exception {
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);

		Range mzRange=new Range(399, 400);
		Range rtInSecRange=new Range(0, 2500000);
		File originalFile=new File("/Volumes/bcsbluessd/kkolotyuk/kkolotyuk_macos_run", "259070_GPF1.dia");
		File newFile=new File("/Volumes/bcsbluessd/kkolotyuk/kkolotyuk_macos_run", "trimmed_259070_GPF1.dia");
		trim(originalFile, newFile, mzRange, rtInSecRange, parameters);
	}
	public static StripeFile trim(File originalFile, File newFile, Range mzRange, Range rtInSecRange, SearchParameters parameters) throws IOException, SQLException, DataFormatException {
		StripeFile stripeFile=new StripeFile();
		stripeFile.openFile();
		
		Logger.logLine("Adding "+originalFile+" to selected from (mz:"+mzRange+" and rt:"+rtInSecRange+")...");
		StripeFileInterface thisStripeFile=StripeFileGenerator.getFile(originalFile, parameters);
		
		HashMap<Range, WindowData> dutyCycleMap=new HashMap<>();
		for (Entry<Range, WindowData> entry : thisStripeFile.getRanges().entrySet()) {
			if (entry.getKey().contains(mzRange.getMiddle())) {
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
			stripes.add(scan);
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
}
