package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import junit.framework.TestCase;

public class GaussianTest extends TestCase {
	public void testGaussian() {
		Distribution g=new Gaussian(0, 1);
		float[] xs=new float[] {-3, -2, -1, 0, 1, 2, 3};
		float[] cdfs=new float[] {0.0013499672813147567f, 0.022750062887256395f, 0.15865526383236372f, 0.500000000f, 0.8413447361676363f, 0.9772499371127437f, 0.9986500327186852f};
		float[] pdfs=new float[] {0.0044318484119380075f, 0.05399096651318806f, 0.24197072451914337f, 0.3989422804014327f, 0.24197072451914337f, 0.05399096651318806f, 0.0044318484119380075f};
		for (int i=0; i<xs.length; i++) {
			assertEquals(cdfs[i], g.getCDF(xs[i]), 0.00001f);
			assertEquals(pdfs[i], g.getPDF(xs[i]), 0.00001f);
		}
	}

}
