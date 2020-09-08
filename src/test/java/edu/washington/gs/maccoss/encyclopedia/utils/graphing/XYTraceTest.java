package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import junit.framework.TestCase;

public class XYTraceTest extends TestCase {
	public void testRound() {
		double[] x = new double[] { 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.49, 2.7, 2.9 };
		double[] y = new double[] { 11, 12, 13, 14, 15, 6, 17, 18, 119, 110, 111, 112, 113, 114, 115 };

		XYTrace trace = new XYTrace(x, y, GraphType.area, "name");
		XYTrace round = XYTrace.round(trace, 1);

		Pair<double[], double[]> xy=round.toArrays();
		float[] expectedx = new float[] { 0.0f, 1.0f, 2.0f, 3.0f };
		float[] expectedy = new float[] { 11.5f, 13.0f, 97.16666666666667f, 114.5f };

		assertEquals(expectedx.length, xy.x.length);
		assertEquals(expectedy.length, xy.y.length);
		for (int i = 0; i < expectedy.length; i++) {
			assertEquals(expectedx[i], xy.x[i], 0.0001);
			assertEquals(expectedy[i], xy.y[i], 0.0001);
		}
	}
}
