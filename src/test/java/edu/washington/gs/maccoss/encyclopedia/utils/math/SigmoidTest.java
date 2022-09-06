package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.Random;

import gnu.trove.list.array.TFloatArrayList;
import junit.framework.TestCase;

public class SigmoidTest extends TestCase {
	public void testSigmoidFitting() {
		TFloatArrayList xList=new TFloatArrayList();
		TFloatArrayList yList=new TFloatArrayList();
		Random random=new Random(7*7*7*7*7);
		
		for (float x = -1; x <= 1; x+=0.01f) {
			xList.add(x);
			yList.add(100f/(1f+(float)Math.exp(-(1+3f*x)))+(float)random.nextGaussian()*10f);
		}
		
		Sigmoid sig=new Sigmoid(1.0f/100000f);
		float[] xs = xList.toArray();
		float[] ys = yList.toArray();
		sig.train(xs, ys, 100000);
		
		assertEquals(0.6340146, sig.getA(), 0.01);
		assertEquals(2.7287354, sig.getB(), 0.1);
		assertEquals(38.410572, sig.getC(), 1);
		assertEquals(57.528805, sig.getD(), 1);
		
		TFloatArrayList xSigList=new TFloatArrayList();
		TFloatArrayList ySigList=new TFloatArrayList();
		for (float x = -2; x <= 2; x+=0.01f) {
			xSigList.add(x);
			ySigList.add(sig.getValue(x));
		}
		
//		XYTrace dataTrace=new XYTrace(xs, ys, GraphType.bighollowpoint, "Data");
//		XYTrace fitTrace=new XYTrace(xSigList.toArray(), ySigList.toArray(), GraphType.line, "Fit");
//		Charter.launchChart("X", "Y", true, dataTrace, fitTrace);
//		
//		try {Thread.sleep(100000);} catch (InterruptedException ie) {};
	}

}
