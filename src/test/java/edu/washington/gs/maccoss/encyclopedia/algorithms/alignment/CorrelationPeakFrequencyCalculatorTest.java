package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import junit.framework.TestCase;

public class CorrelationPeakFrequencyCalculatorTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
	
	public void testPeakFrequencyCalculator() {
		CorrelationPeakFrequencyCalculator calculator=new CorrelationPeakFrequencyCalculator(PARAMETERS.getFragmentTolerance());
		
		calculator.increment(100, 2.0f, 0.1f, true);
		calculator.increment(100.001, 1.0f, 0.1f, true);
		calculator.increment(100.002, 1.0f, 0.1f, true);
		calculator.increment(100.003, 1.0f, 0.1f, true);
		calculator.increment(100.1, 1.0f, 0.1f, true);
		calculator.increment(100.2, 1.0f, 0.1f, true);
		calculator.increment(100.3, 1.0f, 0.1f, true);
		calculator.increment(101, 1.0f, 0.1f, true);
		calculator.increment(101.1, 1.0f, 0.1f, true);
		calculator.increment(101.001, 1.0f, 0.1f, true);
		calculator.increment(101.2, 1.0f, 0.1f, true);
		calculator.increment(101.3, 1.0f, 0.1f, true);
		calculator.increment(501.003, 1.0f, 50f, true);
		calculator.increment(801, 1.0f, 0.1f, true);
		calculator.increment(801.001, 1.0f, 0.1f, true);
		calculator.increment(801.002, 1.0f, 0.1f, true);
		calculator.increment(801.1, 1.0f, 0.1f, true);
		calculator.increment(801.2, 1.0f, 0.1f, true);
		calculator.increment(901.3, 1.0f, 0.1f, true);
		calculator.increment(901.301, 1.0f, 0.1f, true);
		calculator.increment(901.003, 1.0f, 0.1f, true);
		
		double[] masses=calculator.getTopNMasses(4);
		double[] expected=new double[] { 501.003, 801.0, 100.0, 901.3 };
		
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.00001);
		}
		
		//Charter.launchChart("mz", "count", true, new XYTrace(calculator.toPoints(), GraphType.spectrum, "title"));
	}

	public void testPeakFrequencyCalculatorErrors() {
		CorrelationPeakFrequencyCalculator calculator=new CorrelationPeakFrequencyCalculator(PARAMETERS.getFragmentTolerance());
		
		calculator.increment(100, 2.0f, 0.1f, true);
		calculator.increment(100.001, 1.0f, 0.1f, true);
		calculator.increment(100.002, 1.0f, 0.1f, true);
		calculator.increment(100.003, 1.0f, 0.1f, true);
		calculator.increment(100.1, 1.0f, 0.1f, true);
		calculator.increment(100.2, 1.0f, 0.1f, true);
		calculator.increment(100.3, 1.0f, 0.1f, true);
		calculator.increment(101, 1.0f, 0.1f, true);
		calculator.increment(101.1, 1.0f, 0.1f, true);
		calculator.increment(101.001, 1.0f, 0.1f, true);
		calculator.increment(101.2, 1.0f, 0.1f, true);
		calculator.increment(101.3, 1.0f, 0.1f, true);
		
		calculator.increment(601.001, 2.0f, 0.3f, true);
		calculator.increment(601.002, 2.0f, 0.6f, true);
		calculator.increment(601.003, 2.0f, 1f, true);
		calculator.increment(601.003, 0.1f, 0.1f, false);
		calculator.increment(601.003, 0.1f, 0.1f, false);
		calculator.increment(601.003, 0.1f, 0.1f, false);
		
		calculator.increment(501.001, 1.0f, 0.3f, true);
		calculator.increment(501.002, 1.0f, 0.6f, true);
		calculator.increment(501.003, 1.0f, 1f, true);
		calculator.increment(501.003, 0.1f, 0.1f, false);
		calculator.increment(501.003, 0.1f, 0.1f, false);
		calculator.increment(501.003, 0.1f, 0.1f, false);
		calculator.increment(501.003, 0.1f, 0.1f, false);
		
		calculator.increment(801, 1.0f, 0.1f, true);
		calculator.increment(801.001, 1.0f, 0.1f, true);
		calculator.increment(801.002, 1.0f, 0.1f, true);
		calculator.increment(801.1, 1.0f, 0.1f, true);
		calculator.increment(801.2, 1.0f, 0.1f, true);
		calculator.increment(901.3, 1.0f, 0.1f, true);
		calculator.increment(901.301, 1.0f, 0.1f, true);
		calculator.increment(901.003, 1.0f, 0.1f, true);
		
		double[] masses=calculator.getTopNMasses(4);
		double[] expected=new double[] { 601.001, 501.001, 801.0, 100.0 };
		
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.00001);
		}

		calculator.increment(601.003, 10f, 0.1f, false); // one big interference peak kills the transition
		
		masses=calculator.getTopNMasses(4);
		expected=new double[] { 501.001, 801.0, 100.0, 901.3 };

		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.00001);
		}
		
		//Charter.launchChart("mz", "count", true, new XYTrace(calculator.toPoints(), GraphType.spectrum, "title"));
	}
}
