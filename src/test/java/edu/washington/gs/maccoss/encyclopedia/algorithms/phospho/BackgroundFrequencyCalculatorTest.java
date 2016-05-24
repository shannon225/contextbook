package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.MzmlToDIAConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.map.hash.TFloatFloatHashMap;
import junit.framework.TestCase;

public class BackgroundFrequencyCalculatorTest extends TestCase {
	public static void main(String[] args) throws Exception {
		File libraryFile=new File("/Users/searleb/Documents/school/projects/VillenJ_Exactive_HumanPhosphoproteome.elib");
		File diaFile=new File("/Users/searleb/Documents/school/projects/mzml/q06048_rl_MCF7_IMAC_GpX_3.dia");

		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		
		LibraryFile library=new LibraryFile();
		library.openFile(libraryFile);
		
		StripeFileInterface stripefile=MzmlToDIAConverter.getFile(diaFile, parameters);
		
		BackgroundFrequencyCalculator calculator=BackgroundFrequencyCalculator.generateBackground(stripefile, library);

		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();

		ArrayList<String> permutations=PhosphoPermuter.getPermutations("MQS[+80.0]LSLNK", parameters.getAAConstants());
		for (String peptideModSeq : permutations) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
			entryMap.put(peptideModSeq, model);
		}
		
		HashMap<String, FragmentIon[]> uniqueIons=PhosphoLocalizer.getUniqueFragmentIons((byte)2, entryMap, parameters);

		ArrayList<Stripe> stripes=stripefile.getStripes(500.730213, 19.97f*60f-600f, 19.98f*60f+600f, false);
		
		for (Entry<String, FragmentationModel> entry : entryMap.entrySet()) {
			String peptideModSeq=entry.getKey();
			FragmentIon[] targets=uniqueIons.get(peptideModSeq);
			double[] ions=FragmentIon.getMasses(targets);

			float[] frequencies=calculator.getFrequencies(ions, 500.730213, parameters.getFragmentTolerance());
			TFloatFloatHashMap uniqueRtScoreMap=new TFloatFloatHashMap();

			for (Stripe stripe : stripes) {
				float negLogProb=process(parameters, ions, frequencies, stripe);
				System.out.println(peptideModSeq+"\t"+stripe.getScanStartTime()/60f+"\t"+negLogProb);
				uniqueRtScoreMap.put(stripe.getScanStartTime(), negLogProb);
			}

			EValueCalculator uniqueCalculator=new EValueCalculator(uniqueRtScoreMap);
			System.out.println("FINAL: "+peptideModSeq+" --> rt:"+uniqueCalculator.getMaxRT()/60.0f+", s:"+uniqueCalculator.getMaxRawScore()+", e:"+uniqueCalculator.getNegLog10EValue());
		}

		stripes=stripefile.getStripes(500.730213, 19.43f*60f, 19.44f*60f, false);
	}

	private static float process(SearchParameters parameters, double[] ions, float[] frequencies, Stripe stripe) {
		double[] massArray=stripe.getMassArray();
		float logProb=0.0f;
		for (int i=0; i<frequencies.length; i++) {
			float hitProb=frequencies[i]*massArray.length;
			boolean match=parameters.getFragmentTolerance().getIndex(massArray, ions[i]).isPresent();
			if (match) {
				logProb+=Log.log10(hitProb);
			}
		}
		return -logProb;
	}
}
