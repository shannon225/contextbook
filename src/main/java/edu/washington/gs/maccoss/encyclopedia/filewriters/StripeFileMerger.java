package edu.washington.gs.maccoss.encyclopedia.filewriters;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class StripeFileMerger {
	private static final String FILENAME_DELIMITER = ";";
	
	public static void main(String[] args) throws Exception {
		HashMap<String, String> paramMap=PecanParameterParser.getDefaultParameters();
		paramMap.put("-acquisition", "DIA"); // NON-OVERLAPPING!
		SearchParameters parameters=PecanParameterParser.parseParameters(paramMap);
		
		File dir=new File("/Users/searleb/Downloads/wide_iso");
		
		merge(new File[] {new File(dir, "Loo_2020_0406_RJ_102_32.mzML"), new File(dir,"Loo_2020_0406_RJ_103_32.mzML")}, 
				new Range[] {new Range(0, 700), new Range(701,1600)}, 
				new File(dir, "combined.dia"), parameters);
	}

	public static StripeFile merge(File[] fs, Range[] mzRanges, File newFile, SearchParameters parameters) throws IOException, SQLException, DataFormatException {
		return merge(fs, Optional.ofNullable(mzRanges), newFile, parameters);
	}
	public static StripeFile merge(File[] fs, Optional<Range[]> mzRanges, File newFile, SearchParameters parameters) throws IOException, SQLException, DataFormatException {
		StripeFile stripeFile=new StripeFile();
		stripeFile.openFile();
		HashMap<Range, Float> dutyCycleMap=new HashMap<>();
		StringBuilder nameString=new StringBuilder();
		StringBuilder pathString=new StringBuilder();
		
		Arrays.sort(fs);
		
		int scanIndex=0;
		for (int i = 0; i < fs.length; i++) {
			Range mzRange;
			if (mzRanges.isPresent()) {
				mzRange=mzRanges.get()[i];
			} else {
				mzRange=new Range(-Float.MAX_VALUE, Float.MAX_VALUE);
			}
			Logger.logLine("Adding "+fs[i].getName()+" to merged file ("+(i+1)+" of "+fs.length+")...");
			StripeFileInterface thisStripeFile=StripeFileGenerator.getFile(fs[i], parameters);
			dutyCycleMap.putAll(thisStripeFile.getRanges());
			ArrayList<PrecursorScan> precursors = new ArrayList<>();
			for (PrecursorScan scan : thisStripeFile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE)) {
				precursors.add(scan.shallowClone(i, scanIndex++));
			}
			stripeFile.addPrecursor(precursors);
			ArrayList<FragmentScan> stripes = new ArrayList<>();
			for (FragmentScan scan :  thisStripeFile.getStripes(mzRange, -Float.MAX_VALUE, Float.MAX_VALUE, false)) {
				stripes.add(scan.shallowClone(i, scanIndex++));
			}
			stripeFile.addStripe(stripes);
			
			if (nameString.length()>0) nameString.append(FILENAME_DELIMITER);
			nameString.append(fs[i].getName());
			if (pathString.length()>0) pathString.append(FILENAME_DELIMITER);
			pathString.append(fs[i].getAbsolutePath());
			
			thisStripeFile.close();
		}
		Logger.logLine("Finished merging, finalizing "+newFile.getName());

		stripeFile.setFileName(nameString.toString(), null, pathString.toString());
		stripeFile.setRanges(dutyCycleMap);

		stripeFile.saveAsFile(newFile);
		stripeFile.close();
		
		stripeFile=new StripeFile();
		stripeFile.openFile(newFile);
		return stripeFile;
	}
}
