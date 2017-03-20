package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public class BackgroundFrequencyCalculatorTest {
	public static void main(String[] args) throws Exception {
		//File libraryFile=new File("/Users/searleb/Documents/projects/phosphopedia/VillenJ_Exactive_HumanPhosphoproteome.elib");
		//File diaFile=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/dec2015_phospho/110515_bcs_hela_phospho_starved_20mz_500_900.dia");
		//File libraryFile=new File("/Users/searleb/Documents/school/localization_manuscript/hela_phospho/VillenJ_Exactive_HumanPhosphoproteome.elib");
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/hela_phospho/110515_bcs_hela_phospho_starved_20mz_500_900.dia");
		File diaFile=new File("/Users/searleb/Documents/phospho_localization/data/110515_bcs_hela_phospho_starved_20mz_500_900.dia");
		File libraryFile=new File("/Users/searleb/Documents/phospho_localization/data/VillenJ_Exactive_HumanPhosphoproteome.elib");
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		
		//BackgroundFrequencyCalculator calculator=BackgroundFrequencyCalculator.generateBackground(stripefile, library);
		BackgroundFrequencyInterface calculator=BackgroundFrequencyCalculator.generateBackground(stripefile);
		
		Pair<double[], float[]> counterPair=calculator.getRoundedMassCounters(600.0, parameters.getFragmentTolerance());
		double[] masses=counterPair.x;
		float[] counters=counterPair.y;
		for (int i=0; i<counters.length; i++) {
			System.out.println(masses[i]+"\t"+counters[i]);
		}
	}
}
