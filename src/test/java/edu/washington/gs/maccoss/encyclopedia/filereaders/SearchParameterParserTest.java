package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import junit.framework.TestCase;

public class SearchParameterParserTest extends TestCase {
	public void testDefaultParsing() {
		HashMap<String, String> input=PecanParameterParser.getDefaultParameters();
		PecanSearchParameters params=PecanParameterParser.parseParameters(input);

		assertEquals(160.03065030151367, params.getAAConstants().getMass('C'), 0.000001);
		assertEquals(FragmentationType.YONLY, params.getFragType());
		assertEquals(10.0, params.getPrecursorTolerance().getPpmTolerance(), 0.0000001);
		assertEquals(10.0, params.getFragmentTolerance().getPpmTolerance(), 0.0000001);
		assertEquals("Trypsin", params.getEnzyme().getName());
		assertEquals(5, params.getMinPeptideLength());
		assertEquals(100, params.getMaxPeptideLength());
		assertEquals(1, params.getMaxMissedCleavages());
		assertEquals(2, params.getMinCharge());
		assertEquals(3, params.getMaxCharge());
		assertEquals(12, params.getMinEluteTime());
		assertEquals(1, params.getNumberOfReportedPeaks());
		assertEquals(false, params.isAddDecoysToBackgound());
	}

	public void testParsing() {
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("-fixed", "C=3,T=6");
		map.put("-frag", "CID");
		map.put("-ptol", "12");
		map.put("-ftol", "13");
		map.put("-enzyme", "Lys-C");
		map.put("-minLength", "6");
		map.put("-maxLength", "40");
		map.put("-maxMissedCleavage", "2");
		map.put("-minCharge", "1");
		map.put("-maxCharge", "4");
		map.put("-minEluteTime", "20");
		map.put("-numberOfReportedPeaks", "3");
		map.put("-addDecoysToBackground", "true");
		PecanSearchParameters params=PecanParameterParser.parseParameters(map);
		
		assertEquals(103.009185+3.0, params.getAAConstants().getMass('C'), 0.000001);
		assertEquals(101.047679+6.0, params.getAAConstants().getMass('T'), 0.000001);
		assertEquals(FragmentationType.CID, params.getFragType());
		assertEquals(12.0, params.getPrecursorTolerance().getPpmTolerance(), 0.0000001);
		assertEquals(13.0, params.getFragmentTolerance().getPpmTolerance(), 0.0000001);
		assertEquals("Lys-C", params.getEnzyme().getName());
		assertEquals(6, params.getMinPeptideLength());
		assertEquals(40, params.getMaxPeptideLength());
		assertEquals(2, params.getMaxMissedCleavages());
		assertEquals(1, params.getMinCharge());
		assertEquals(4, params.getMaxCharge());
		assertEquals(20, params.getMinEluteTime());
		assertEquals(3, params.getNumberOfReportedPeaks());
		assertEquals(true, params.isAddDecoysToBackgound());
	}

}
