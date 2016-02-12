package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import com.google.common.base.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class EncyclopediaToBLIBTest {
	public static void main(String[] args) throws Exception {
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-deconvoluteOverlappingWindows", "true");
		map.put("-fixed", "");
		SearchParameters parameters=SearchParameterParser.parseParameters(map);
		
		boolean generateQuantLib=true;
		File blibFile;
		SearchJobData job1;
		LibraryFile libraryTemplate=null;
		if (generateQuantLib) {
			blibFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000_Quant.blib");
			job1=getData(parameters, "/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.mzML");

			File libraryFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.elib");
			libraryTemplate=new LibraryFile();
			libraryTemplate.openFile(libraryFile);
		} else {
			blibFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.blib");
			job1=getData(parameters, "/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.mzML");
		}
		
		//File blibFile=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz.blib");
		//SearchJobData job1=getData(parameters, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110515_bcs_hela_phospho_starved_20mz_500_900.mzML");
		
		/*
		SearchJobData job1=getData(parameters, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_500_660.mzML");
		SearchJobData job2=getData(parameters, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_660_820.mzML");
		SearchJobData job3=getData(parameters, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_820_980.mzML");
		SearchJobData job4=getData(parameters, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_980_1140.mzML");
		SearchJobData job5=getData(parameters, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_1140_1300.mzML");
		 */
		ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
		jobs.add(job1);
		/*
		jobs.add(job2);
		jobs.add(job3);
		jobs.add(job4);
		jobs.add(job5);
		*/
		
		SearchToBLIB.convert(new EmptyProgressIndicator(), jobs, blibFile, Optional.fromNullable(libraryTemplate));
	}

	private static SearchJobData getData(SearchParameters parameters, String dia) {
		File diaFile=new File(dia);
		File outputFile=new File(diaFile.getAbsolutePath()+".percolator.txt");
		File featureFile=new File(diaFile.getAbsolutePath()+".features.txt");
		SearchJobData job=new SearchJobData(diaFile, featureFile, outputFile, parameters, "1");
		return job;
	}
}
