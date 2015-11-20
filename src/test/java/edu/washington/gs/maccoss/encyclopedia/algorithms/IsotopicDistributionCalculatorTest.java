package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import junit.framework.TestCase;

public class IsotopicDistributionCalculatorTest extends TestCase {
	private static final SearchParameters PARAMETERS=new SearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));

	public void testGetIsotopeDistributionLongString() {
		float[] dist=IsotopicDistributionCalculator.getIsotopeDistribution("LACDEFQFEDCAIRLACDEFQFEDCAIR", PARAMETERS.getAAConstants());
		float[] expected=new float[] {0.49305102f, 0.92382437f, 1.0f, 0.7919198f, 0.5046976f};
		for (int i=0; i<dist.length; i++) {
			assertEquals(expected[i], dist[i], 0.00001f);
		}
	}

	public void testGetIsotopeDistributionString() {
		float[] dist=IsotopicDistributionCalculator.getIsotopeDistribution("LACDEFQFEDCAIR", PARAMETERS.getAAConstants());
		float[] expected=new float[] {1.0f, 0.9371813f, 0.576593f, 0.26531103f, 0.09955499f};
		for (int i=0; i<dist.length; i++) {
			assertEquals(expected[i], dist[i], 0.00001f);
		}
	}

	public void testGetIsotopeDistributionIntArray() {
		int[] freq=new int[] {20, 20, 4, 8, 2};
		float[] dist=IsotopicDistributionCalculator.getIsotopeDistribution(freq);
		float[] expected=new float[] {1.0f, 0.26866445f, 0.13181677f, 0.028259907f, 0.0063654757f};
		for (int i=0; i<dist.length; i++) {
			assertEquals(expected[i], dist[i], 0.00001f);
		}
	}

}
