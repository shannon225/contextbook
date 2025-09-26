package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import junit.framework.TestCase;

public class LinearInterpolatedFunctionTest extends TestCase {
	public void testInterpolation() {
		double[] x = {
			    -1.571, -1.471, -1.371, -1.271, -1.171, -1.071,
			    -0.971, -0.871, -0.771, -0.671, -0.571, -0.471,
			    -0.371, -0.271, -0.171, -0.071,  0.029,  0.129,
			     0.229,  0.329,  0.429,  0.529,  0.629,  0.729,
			     0.829,  0.929,  1.029,  1.129,  1.229,  1.329,
			     1.429,  1.529
			};

		double[] y = {
			    0.000, 0.002, 0.010, 0.022, 0.039, 0.061, 0.087, 0.118,
			    0.152, 0.189, 0.230, 0.273, 0.319, 0.366, 0.415, 0.465,
			    0.515, 0.564, 0.614, 0.662, 0.708, 0.752, 0.794, 0.833,
			    0.869, 0.901, 0.928, 0.952, 0.971, 0.985, 0.995, 1.000
			};
		
		ArrayList<XYPoint> points=XYTrace.toPoints(x, y);
		
		float[] expectedX = {
			    -2.0f, -1.8f, -1.6f, -1.4f, -1.2f, -1.0f,
			    -0.8f, -0.6f, -0.4f, -0.2f,  0.0f,  0.2f,
			     0.4f,  0.6f,  0.8f,  1.0f,  1.2f,  1.4f,
			     1.6f,  1.8f,  2.0f
			};
		
		float[] expectedY = {
			    0.0f, 0.0f, 0.000f, 0.007f, 0.034f, 0.079f,
			    0.141f, 0.218f, 0.305f, 0.401f, 0.500f, 0.599f,
			    0.695f, 0.782f, 0.859f, 0.921f, 0.966f, 0.993f,
			    1.000f, 1.000f, 1.000f
			};
		
		LinearInterpolatedFunction function=new LinearInterpolatedFunction(points);
		for (int i = 0; i < expectedX.length; i++) {
			assertEquals(expectedY[i], function.getYValue(expectedX[i]), 0.002f);
		}
	}
}
