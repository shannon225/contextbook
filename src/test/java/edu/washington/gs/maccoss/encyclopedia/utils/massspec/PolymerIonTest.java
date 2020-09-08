package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import junit.framework.TestCase;

public class PolymerIonTest extends TestCase {

	public void testPEG() {
		double[] expected = new double[] { 63.044055898362096, 107.07027064620769, 151.0964853940533, 195.1227001418989,
				239.14891488974447, 283.17512963759015, 327.20134438543573, 371.2275591332813, 415.2537738811269,
				459.2799886289725, 503.30620337681813, 547.3324181246637, 591.3586328725092, 635.3848476203549,
				679.4110623682004, 723.4372771160461, 767.4634918638917, 811.4897066117372, 855.5159213595829,
				899.5421361074284, 943.568350855274, 987.5945656031197 };

		ArrayList<PolymerIon> ions = PolymerIon.getPolymerProducts(Polymer.peg, new Range(0f, 1000f));
		
		assertEquals(expected.length, ions.size());
		for (int i = 0; i < expected.length; i++) {
			assertEquals(expected[i], ions.get(i).getMass(), 0.00001);
			assertEquals(i+1, ions.get(i).getN());
		}
	}
}
