package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import junit.framework.TestCase;

public class CorrelationPeakFrequencyCalculatorTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false);
	
	public void testPeakFrequencyCalculator() {
		CorrelationPeakFrequencyCalculator calculator=new CorrelationPeakFrequencyCalculator(PARAMETERS.getFragmentTolerance());
		
		calculator.increment(100, 2.0f, 0.1f);
		calculator.increment(100.001, 1.0f, 0.1f);
		calculator.increment(100.002, 1.0f, 0.1f);
		calculator.increment(100.003, 1.0f, 0.1f);
		calculator.increment(100.1, 1.0f, 0.1f);
		calculator.increment(100.2, 1.0f, 0.1f);
		calculator.increment(100.3, 1.0f, 0.1f);
		calculator.increment(101, 1.0f, 0.1f);
		calculator.increment(101.1, 1.0f, 0.1f);
		calculator.increment(101.001, 1.0f, 0.1f);
		calculator.increment(101.2, 1.0f, 0.1f);
		calculator.increment(101.3, 1.0f, 0.1f);
		calculator.increment(501.003, 1.0f, 50f);
		calculator.increment(801, 1.0f, 0.1f);
		calculator.increment(801.001, 1.0f, 0.1f);
		calculator.increment(801.002, 1.0f, 0.1f);
		calculator.increment(801.1, 1.0f, 0.1f);
		calculator.increment(801.2, 1.0f, 0.1f);
		calculator.increment(901.3, 1.0f, 0.1f);
		calculator.increment(901.301, 1.0f, 0.1f);
		calculator.increment(901.003, 1.0f, 0.1f);
		
		double[] masses=calculator.getTopNMasses(4);
		double[] expected=new double[] { 501.003, 801.0, 100.0, 901.3 };
		
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.00001);
		}
		
		//Charter.launchChart("mz", "count", true, new XYTrace(calculator.toPoints(), GraphType.spectrum, "title"));
	}

}
