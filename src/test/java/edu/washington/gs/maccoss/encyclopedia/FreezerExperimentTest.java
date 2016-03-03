package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class FreezerExperimentTest {
	public static void main(String[] args) throws Exception {
		// process all files in this directory
		File dir=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/freezer/");
		
		HashMap<String, String> map1=SearchParameterParser.getDefaultParameters();
		map1.put("-deconvoluteOverlappingWindows", "true");
		SearchParameters searchParameters=SearchParameterParser.parseParameters(map1);

		HashMap<String, String> map2=SearchParameterParser.getDefaultParameters();
		map2.put("-deconvoluteOverlappingWindows", "true");
		map2.put("-fixed", "");
		SearchParameters extractionParameters=SearchParameterParser.parseParameters(map2);
		
		EncyclopediaOneScoringFactory factory=new EncyclopediaOneScoringFactory(searchParameters);
		File libraryFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/hela_6mz.elib");

		File libraryTemplateFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.elib");
		LibraryFile libraryTemplate=new LibraryFile();
		libraryTemplate.openFile(libraryTemplateFile);
		
		File[] files=dir.listFiles();
		files=new File[] {
				new File("/Users/searleb/Documents/projects/encyclopedia/mzml/freezer/121115_bcs_hela_24mz_400_1000_rt_7D_1.mzML")
		};
		long totalTime=0;
		int fileCount=0;
		
		int totalCount=0;
		for (File file : files) {
			if (file.getName().endsWith("mzML")) {
				totalCount++;
			}
		}
		for (File file : files) {
			if (file.getName().endsWith("mzML")) {
				long currentTime=System.currentTimeMillis();
				EncyclopediaTest.run(file, libraryFile, factory);

				File blibFile=new File(file.getAbsolutePath()+".quant.blib");
				SearchJobData job=getData(extractionParameters, file);

				ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
				jobs.add(job);

				SearchToBLIB.convert(new EmptyProgressIndicator(), jobs, blibFile, Optional.ofNullable(libraryTemplate));
				
				fileCount++;
				long fileTime=System.currentTimeMillis()-currentTime;
				totalTime+=fileTime;
				float averageTimePer=totalTime/1000f/60f/fileCount;
				float remaining=(totalCount-fileCount)*averageTimePer/60f;
				System.out.println("Processed "+fileCount+"/"+totalCount+" files in "+Math.round(10f*totalTime/1000f/60f/60f)/10f+" hours (average of "+Math.round(10f*averageTimePer)/10f+" minutes per file, "+Math.round(10f*remaining)/10f+" hours remaining)");
			}
		}
	}
	
	private static SearchJobData getData(SearchParameters parameters, File diaFile) {
		File outputFile=new File(diaFile.getAbsolutePath()+".percolator.txt");
		File featureFile=new File(diaFile.getAbsolutePath()+".features.txt");
		SearchJobData job=new SearchJobData(diaFile, featureFile, outputFile, parameters, "1");
		return job;
	}
}
