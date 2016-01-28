package edu.washington.gs.maccoss.encyclopedia;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;

public class EncyclopediaTest {
	public static void main(String[] args) {
		/*
		Encyclopedia.main(new String[] { "-i", "/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.mzML", 
				"-l", "/Users/searleb/Documents/projects/encyclopedia/mzml/cptac2_human_hcd_selected.elib",
				//"-l", "/Users/searleb/Documents/projects/encyclopedia/mzml/hela_6mz.elib",
				"-targetWindowCenter", "750", 
				"-deconvoluteOverlappingWindows", "true"});
		*/

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
				
	}
}
