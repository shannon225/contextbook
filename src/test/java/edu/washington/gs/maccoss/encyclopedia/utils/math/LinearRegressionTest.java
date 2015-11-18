package edu.washington.gs.maccoss.encyclopedia.utils.math;

import junit.framework.TestCase;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public class LinearRegressionTest extends TestCase {
	/**
	 * test example from: https://www.clemson.edu/ces/phoenix/tutorials/excel/regression.html
	 */
	public void testRegression() {
		float[] x=new float[] {1f, 2.3f, 3.1f, 4.8f, 5.6f, 6.3f};
		float[] y=new float[] {2.6f, 2.8f, 3.1f, 4.7f, 5.1f, 5.3f};
		Pair<Float, Float> equation=LinearRegression.getRegression(x, y);
		assertEquals(0.5842f, equation.x, 0.0001f);
		assertEquals(1.6842f, equation.y, 0.0001f);
	}
}
