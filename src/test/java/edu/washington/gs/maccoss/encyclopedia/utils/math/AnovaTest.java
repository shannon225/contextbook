package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;

import org.apache.commons.math3.stat.inference.TestUtils;

import junit.framework.TestCase;

public class AnovaTest extends TestCase {
	public void testAnova() {
		double[] a = { 0, 0, 4478.8267, 0, 0, 3598.408 };
		double[] b = { 0, 0, 0, 0, 0, 0 };
		double[] c = { 0, 0, 0, 0, 0, 0 };
		double[] d = { 16935.947, 17843.93, 16340.3125, 11211.446, 22069.186, 14507.263 };
		double[] e = { 17486.436, 18743.143, 15767.936, 24486.268, 14517.888, 18481.402 };
		double[] f = { 2164.6895, 0, 0, 889.359, 0, 0 };
		
		double[][] data=new double[][] {a, b, c, d, e, f};

		ArrayList<double[]> classes=new ArrayList<>();
		for (int i=0; i<data.length; i++) {
			classes.add(data[i]);
		}
		System.out.println(TestUtils.oneWayAnovaPValue(classes));
	}
}
