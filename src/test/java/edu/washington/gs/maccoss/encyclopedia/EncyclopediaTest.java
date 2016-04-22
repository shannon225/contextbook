package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class EncyclopediaTest {
	public static void main(String[] args) throws Exception {
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-deconvoluteOverlappingWindows", "true");
		//map.put("-targetWindowCenter", "750");
		SearchParameters parameters=SearchParameterParser.parseParameters(map);
		EncyclopediaOneScoringFactory factory=new EncyclopediaOneScoringFactory(parameters);
		
		/*
		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/encyclopedia/mzml/momo/20150301_A1_DIA_1.mzML", 
				"-l", "/Users/searleb/Documents/projects/encyclopedia/mzml/momo/20151125_122sample_id_iRT.elib",
				"-deconvoluteOverlappingWindows", "true"});
		*/

		File libraryFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/hela_6mz.elib");
		run(new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.mzML"), libraryFile, factory);
		
		/*
		libraryFile=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib");

	
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_500_660.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_660_820.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_820_980.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_980_1140.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_1140_1300.mzML"), libraryFile, factory);
		
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_starved_6mz_500_660.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_starved_6mz_660_820.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_starved_6mz_820_980.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_starved_6mz_980_1140.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_starved_6mz_1140_1300.mzML"), libraryFile, factory);

		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110515_bcs_hela_phospho_starved_20mz_500_900.mzML"), libraryFile, factory);
		run(new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110515_bcs_hela_phospho_starved_20Mz_900_1300.mzML"), libraryFile, factory);
		*/
		
		/*
		// NEGATIVE!
		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/encyclopedia/mzml/yeast/Q_2014_0523_12_0_amol_uL_20mz.mzML", 
				//"-l", "/Users/searleb/Documents/projects/encyclopedia/mzml/cptac2_human_hcd_selected.elib",
				"-l", "/Users/searleb/Documents/projects/encyclopedia/mzml/hela_6mz.elib",
				//"-targetWindowCenter", "750", 
				"-deconvoluteOverlappingWindows", "true"});
		
		/*
		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_980_1140.mzML", 
				"-l", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib",
				"-deconvoluteOverlappingWindows", "true"});
		
		/*
		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_500_660.mzML", 
				"-l", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib",
				"-deconvoluteOverlappingWindows", "true"});

		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_660_820.mzML", 
				"-l", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib",
				"-deconvoluteOverlappingWindows", "true"});

		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_820_980.mzML",
				"-l", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib",
				"-deconvoluteOverlappingWindows", "true"});

		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_980_1140.mzML", 
				"-l", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib",
				"-deconvoluteOverlappingWindows", "true"});

		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_1140_1300.mzML", 
				"-l", "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/VillenJ_Exactive_HumanPhosphoproteome.elib",
				"-deconvoluteOverlappingWindows", "true"});
		
		HashMap<String, String> defaultParameters=PecanParameterParser.getDefaultParameters();
		defaultParameters.put("-frag", "CID");
		PecanSearchParameters parameters=PecanParameterParser.parseParameters(defaultParameters);
		File fastaFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/UP000005640_9606.fasta");
		File blibFile=new File("/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/hela_phospho_6mz.blib");
		ArrayList<FastaEntry> targets=null;
		
		PecanJobData job1=PecanToBLIB.getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_500_660.mzML");
		PecanJobData job2=PecanToBLIB.getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_660_820.mzML");
		PecanJobData job3=PecanToBLIB.getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_820_980.mzML");
		PecanJobData job4=PecanToBLIB.getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_980_1140.mzML");
		PecanJobData job5=PecanToBLIB.getData(parameters, fastaFile, targets, "/Users/searleb/Documents/projects/pecan/bcs_hela/phospho/110415_bcs_hela_phospho_igf1_6mz_1140_1300.mzML");

		ArrayList<PecanJobData> jobs=new ArrayList<PecanJobData>();
		jobs.add(job1);
		jobs.add(job2);
		jobs.add(job3);
		jobs.add(job4);
		jobs.add(job5);
		PecanToBLIB.convert(new EmptyProgressIndicator(), jobs, blibFile);
		*/	
	}

	public static void run(File diaFile, File libraryFile, EncyclopediaOneScoringFactory factory) throws IOException, SQLException, DataFormatException, ExecutionException, InterruptedException {
		LibraryInterface library=BlibToLibraryConverter.getFile(libraryFile);
		EncyclopediaJobData job;
		job=new EncyclopediaJobData(diaFile, library, factory);
		Encyclopedia.runSearch(new EmptyProgressIndicator(), job);
	}
}
