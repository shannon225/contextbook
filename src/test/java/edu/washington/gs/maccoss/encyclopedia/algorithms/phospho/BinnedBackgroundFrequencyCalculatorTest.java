package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileGenerator;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import junit.framework.TestCase;

public class BinnedBackgroundFrequencyCalculatorTest extends TestCase {
	public static void main(String[] args) throws Exception {
		//File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/hela_phospho/110515_bcs_hela_phospho_starved_20mz_500_900.dia");
		File diaFile=new File("/Users/searleb/Documents/school/localization_manuscript/elibs/mcf7/22jun2016_mcf7_phospho_1a.dia");
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		StripeFileInterface stripefile=StripeFileGenerator.getFile(diaFile, parameters);
		
		BackgroundFrequencyInterface calculator=BinnedBackgroundFrequencyCalculator.generateBackground(stripefile);
		
		Pair<double[], float[]> counterPair=calculator.getRoundedMassCounters(600.0, parameters.getFragmentTolerance());
		double[] masses=counterPair.x;
		float[] counters=counterPair.y;
		for (int i=0; i<counters.length; i++) {
			if (counters[i]>0.0f) {
				System.out.println(masses[i]+"\t"+counters[i]);
			}
		}
	}
	
	public void testIndexing() {
		for (int i=0; i<1000000; i++) {
			double mass=Math.random()*1500.0;
			int index=BinnedBackgroundFrequencyCalculator.getIndex(mass);
			double retrievedMass=BinnedBackgroundFrequencyCalculator.getMass(index);
			assertEquals(mass, retrievedMass, 0.1);
			int retrieveIndex=BinnedBackgroundFrequencyCalculator.getIndex(mass);
			assertEquals(index, retrieveIndex);
		}
	}

}
