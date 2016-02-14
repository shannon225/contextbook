package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;

public class PrecursorIntegrator {
	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();

		File dir=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/freezer/");
		File[] files=dir.listFiles();
		for (File file : files) {
			if (file.getName().endsWith("mzML")) {
				float tic=0.0f;
				StripeFileInterface stripefile=MzmlToDIAConverter.getFile(file, parameters);
				ArrayList<PrecursorScan> scans=stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE);
				for (PrecursorScan scan : scans) {
					float[] intensities=scan.getIntensityArray();
					for (int i=0; i<intensities.length; i++) {
						tic+=intensities[i];
					}
				}
				System.out.println(file.getName()+"\t"+tic);
			}
		}
	}
}
