package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import junit.framework.TestCase;

public class PhosphoPermuterTest extends TestCase {
	public void testPermutations() {
		// normal test one site
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<String> permutations=PhosphoPermuter.getPermutations("AGS[+80]VANAR", PeptideModification.phosphorylation, parameters.getAAConstants());
		String[] expected=new String[] { "AGS[+79.96633]VANAR" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}
		
		// normal test
		parameters=SearchParameterParser.getDefaultParametersObject();
		permutations=PhosphoPermuter.getPermutations("SGS[+80]VSNYR", PeptideModification.phosphorylation, parameters.getAAConstants());
		expected=new String[] { "S[+79.96633]GSVSNYR", "SGS[+79.96633]VSNYR", "SGSVS[+79.96633]NYR", "SGSVSNY[+79.96633]R" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}
		
		// multiple mods
		permutations=PhosphoPermuter.getPermutations("AHT[+80]FSHPPS[+80]STKR", PeptideModification.phosphorylation, parameters.getAAConstants());
		expected=new String[] { "AHT[+79.96633]FS[+79.96633]HPPSSTKR", "AHT[+79.96633]FSHPPS[+79.96633]STKR", "AHT[+79.96633]FSHPPSS[+79.96633]TKR", "AHT[+79.96633]FSHPPSST[+79.96633]KR",
				"AHTFS[+79.96633]HPPS[+79.96633]STKR", "AHTFS[+79.96633]HPPSS[+79.96633]TKR", "AHTFS[+79.96633]HPPSST[+79.96633]KR", "AHTFSHPPS[+79.96633]S[+79.96633]TKR",
				"AHTFSHPPS[+79.96633]ST[+79.96633]KR", "AHTFSHPPSS[+79.96633]T[+79.96633]KR" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}

		// fully modified
		permutations=PhosphoPermuter.getPermutations("S[+80]GS[+80]VS[+80]NY[+80]R", PeptideModification.phosphorylation, parameters.getAAConstants());
		expected=new String[] { "S[+79.96633]GS[+79.96633]VS[+79.96633]NY[+79.96633]R" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}
		
		// doesn't affect modifications
		permutations=PhosphoPermuter.getPermutations("C[+57]SS[+80]VTGVQR", PeptideModification.phosphorylation, parameters.getAAConstants());
		expected=new String[] { "C[+57.0214635]S[+79.96633]SVTGVQR", "C[+57.0214635]SS[+79.96633]VTGVQR", "C[+57.0214635]SSVT[+79.96633]GVQR" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}
	}
	
	public void testFragmentation() {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<String> permutations=PhosphoPermuter.getPermutations("SGS[+80]VSNYR", PeptideModification.phosphorylation, parameters.getAAConstants());
		double[][] expected=new double[][] {new double[] {
				// S[+79.96633]GSVSNYR
				175.11895217, //y1
				338.18227217000003, //y2
				452.22519917000005, //y3
				539.2572271700001, //y4
				638.32564117, //y5
				725.35766917, //y6
				782.37913317, //y7
				851.40059517, //y8-NL
				949.37749117, //y8
			}, new double[] {
				// SGS[+79.96633]VSNYR
				175.11895217, //y1
				338.18227217000003, //y2
				452.22519917000005, //y3
				539.2572271700001, //y4
				638.32564117, //y5
				707.34710317, //y6-NL
				764.36856717, //y7-NL
				805.32399917, //y6
				851.40059517, //y8-NL
				862.34546317, //y7
				949.37749117, //y8
			}, new double[] {
				// SGSVS[+79.96633]NYR
				175.11895217, //y1
				338.18227217000003, //y2
				452.22519917000005, //y3
				521.24666117, //y4-NL
				619.22355717, //y4
				620.31507517, //y5-NL
				707.34710317, //y6-NL
				718.29197117, //y5
				764.36856717, //y7-NL
				805.32399917, //y6
				851.40059517, //y8-NL
				862.34546317, //y7
				949.37749117, //y8
			}, new double[] {
				// SGSVSNY[+79.96633]R // no neutral loss
				175.11895217, //y1
				418.14860217, //y2
				532.19152917, //y3
				619.2235571699999, //y4
				718.2919711699999, //y5
				805.3239991699999, //y6
				862.3454631699999, //y7
				949.3774911699999, //y8
		}};

		for (int i=0; i<permutations.size(); i++) {
			FragmentationModel model=new FragmentationModel(permutations.get(i), parameters.getAAConstants());
			FragmentIon[] yions=model.getYIons();
			for (int j=0; j<yions.length; j++) {
				assertEquals(expected[i][j], yions[j].mass, 0.01);
			}
		}
	}
}
