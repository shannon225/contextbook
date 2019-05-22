package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import junit.framework.TestCase;

public class SearchParametersTest extends TestCase {
	public void testReadParameters() throws Exception {
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-runPhosphoLocalization", "true");
		map.put("-deconvoluteOverlappingWindows", "true");
		map.put("-poffset", "1000"); // definitely not default!
		
		SearchParameters params=SearchParameterParser.parseParameters(map);
		assertEquals(1000.0, params.getPrecursorOffsetPPM());
		params.savePreferences(null, null);
		
		SearchParameters readParams=SearchParameterParser.parseParameters(SearchParameters.readPreferences());
		assertEquals(1000.0, readParams.getPrecursorOffsetPPM());
		
		map.put("-poffset", "-1000"); // definitely not default!
		
		params=SearchParameterParser.parseParameters(map);
		assertEquals(-1000.0, params.getPrecursorOffsetPPM());
		params.savePreferences(null, null);
		
		readParams=SearchParameterParser.parseParameters(SearchParameters.readPreferences());
		assertEquals(-1000.0, readParams.getPrecursorOffsetPPM());
	}

	public void testEntrapmentFDR() {
		float FDR=0.01f;
		for (float i = 0.0f; i <= 2.05f; i+=0.1f) {
			System.out.println(i+"\t"+SearchParameters.getEffectivePercolatorThreshold(0.01f, i)+"\t"+(FDR * (1+i)*(1-(((1+i)-1)*FDR)))+"\t"+(FDR*(1+i)));
		}
	}
}
