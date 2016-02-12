package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import com.google.common.base.Optional;

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
		for (File file : files) {
			if (file.getName().endsWith("mzML")) {
				EncyclopediaTest.run(file, libraryFile, factory);

				File blibFile=new File(file.getAbsolutePath()+".quant.blib");
				SearchJobData job=getData(extractionParameters, file);

				ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
				jobs.add(job);

				SearchToBLIB.convert(new EmptyProgressIndicator(), jobs, blibFile, Optional.fromNullable(libraryTemplate));
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
