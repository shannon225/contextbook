package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import junit.framework.TestCase;

public class PhosphoPermuterTest extends TestCase {
	public void testPermutations() {
		// normal test
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<String> permutations=PhosphoPermuter.getPermutations("SGS[+80]VSNYR", parameters.getAAConstants());
		String[] expected=new String[] { "S[+79.96633]GSVSNYR", "SGS[+79.96633]VSNYR", "SGSVS[+79.96633]NYR", "SGSVSNY[+79.96633]R" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}
		
		// multiple mods
		permutations=PhosphoPermuter.getPermutations("AHT[+80]FSHPPS[+80]STKR", parameters.getAAConstants());
		expected=new String[] { "AHT[+79.96633]FS[+79.96633]HPPSSTKR", "AHT[+79.96633]FSHPPS[+79.96633]STKR", "AHT[+79.96633]FSHPPSS[+79.96633]TKR", "AHT[+79.96633]FSHPPSST[+79.96633]KR",
				"AHTFS[+79.96633]HPPS[+79.96633]STKR", "AHTFS[+79.96633]HPPSS[+79.96633]TKR", "AHTFS[+79.96633]HPPSST[+79.96633]KR", "AHTFSHPPS[+79.96633]S[+79.96633]TKR",
				"AHTFSHPPS[+79.96633]ST[+79.96633]KR", "AHTFSHPPSS[+79.96633]T[+79.96633]KR" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}

		// fully modified
		permutations=PhosphoPermuter.getPermutations("S[+80]GS[+80]VS[+80]NY[+80]R", parameters.getAAConstants());
		expected=new String[] { "S[+79.96633]GS[+79.96633]VS[+79.96633]NY[+79.96633]R" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}
		
		// doesn't affect modifications
		permutations=PhosphoPermuter.getPermutations("C[+57]SS[+80]VTGVQR", parameters.getAAConstants());
		expected=new String[] { "C[+57.0]S[+79.96633]SVTGVQR", "C[+57.0]SS[+79.96633]VTGVQR", "C[+57.0]SSVT[+79.96633]GVQR" };
		assertEquals(expected.length, permutations.size());
		for (int i=0; i<permutations.size(); i++) {
			assertEquals(expected[i], permutations.get(i));
		}
	}
}
