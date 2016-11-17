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
		params.savePreferences();
		
		SearchParameters readParams=SearchParameters.readPreferences();
		assertEquals(1000.0, readParams.getPrecursorOffsetPPM());
		
		map.put("-poffset", "-1000"); // definitely not default!
		
		params=SearchParameterParser.parseParameters(map);
		assertEquals(-1000.0, params.getPrecursorOffsetPPM());
		params.savePreferences();
		
		readParams=SearchParameters.readPreferences();
		assertEquals(-1000.0, readParams.getPrecursorOffsetPPM());
	}

}
