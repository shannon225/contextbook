package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class PecanToBLIBTest {
	public static void main(String[] args) {
		long time=System.currentTimeMillis();
		
		HashMap<String, String> defaultParameters=PecanParameterParser.getDefaultParameters();
		defaultParameters.put("-frag", "CID");
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(defaultParameters);
		File fastaFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/UP000005640_9606.fasta");
		File blibFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/slow_test_hela_6mz.blib");
		ArrayList<FastaEntry> targets=null;
		
		/*SearchJobData job1=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/encyclopedia/mzml/121015_BCS_HeLa_6mz_400_500.mzML");
		SearchJobData job2=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/encyclopedia/mzml/121015_BCS_HeLa_6mz_500_600.mzML");
		SearchJobData job3=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/encyclopedia/mzml/121015_BCS_HeLa_6mz_600_700.mzML");
		SearchJobData job4=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/encyclopedia/mzml/121015_BCS_HeLa_6mz_700_800.mzML");
		SearchJobData job5=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/encyclopedia/mzml/121015_BCS_HeLa_6mz_800_900.mzML");
		SearchJobData job6=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/encyclopedia/mzml/121015_BCS_HeLa_6mz_900_1000.mzML");

		ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
		jobs.add(job1);
		jobs.add(job2);
		jobs.add(job3);
		jobs.add(job4);
		jobs.add(job5);
		jobs.add(job6);*/
		
		SearchJobData job1=getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/encyclopedia/mzml/freezer/121115_bcs_hela_24mz_400_1000_0D_1.mzML");

		ArrayList<SearchJobData> jobs=new ArrayList<SearchJobData>();
		jobs.add(job1);
		
		LibraryInterface libraryTemplate=null;
		SearchToBLIB.convert(new EmptyProgressIndicator(), jobs, blibFile, true, Optional.ofNullable(libraryTemplate));
		
		System.out.println((System.currentTimeMillis()-time)/1000+" seconds");
	}

	private static SearchJobData getData(PecanSearchParameters parameters, File fastaFile, ArrayList<FastaEntry> targets, String dia) {
		File diaFile1=new File(dia);
		File outputFile1=new File(diaFile1.getAbsolutePath()+".percolator.txt");
		File featureFile1=new File(outputFile1.getAbsolutePath()+".features.txt");
		PecanScoringFactory factory1=new PecanOneScoringFactory(parameters, featureFile1);
		SearchJobData job1=new PecanJobData(Optional.ofNullable(targets), diaFile1, fastaFile,featureFile1, outputFile1, factory1);
		return job1;
	}
}
