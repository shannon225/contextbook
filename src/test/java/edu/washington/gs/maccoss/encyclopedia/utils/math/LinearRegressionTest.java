package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import junit.framework.TestCase;

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
	
	public void testBoundary() {
		LinearRegression r=new LinearRegression(new float[] {0, 1}, new float[] {0, 1});
		assertEquals(0.0f, r.b);
		assertEquals(1.0f, r.m);
		r=new LinearRegression(new float[] {1}, new float[] {1});
		assertEquals(1f, r.b);
		assertEquals(0f, r.m);
	}
	
	public void testBoundary2() {
		LinearRegression r=new LinearRegression(new float[] {0, 1}, new float[] {0, 1}, new XYPoint(0f, 0.1f));
		assertEquals(0.1f, r.b);
		assertEquals(0.9f, r.m);

		r=new LinearRegression(new float[] {0, 1}, new float[] {0, 1}, new XYPoint(0.1f, 0.1f));
		assertEquals(0f, r.b);
		assertEquals(1f, r.m);
	}
	
	public void testRegressionWithFixedIntercept() {
		float[] x = { -0.6661775f, -0.3309932f, -0.16694789f, 0.0f };
		float[] y = { -0.98749655f, -0.5642714f, -0.2620367f, 0.0f };
		
		Pair<Float, Float> r=LinearRegression.getRegressionWithFixedIntercept(x, y, new XYPoint(-1.0, -2.0628977));
		assertEquals(2.1858542f, r.x, 000001f);
		assertEquals(0.1229565f, r.y, 000001f);
	}
}
