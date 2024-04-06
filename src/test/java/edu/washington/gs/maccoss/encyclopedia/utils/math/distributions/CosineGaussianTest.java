package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import junit.framework.TestCase;

public class CosineGaussianTest extends TestCase {
	public static void main(String[] args) {
		CosineGaussian dist=new CosineGaussian(0, 1, 1);
		Gaussian dist2=new Gaussian(0, 1, 1);
		Range range=new Range(-3, 3);
		for (int i = 0; i < 10; i++) {
			System.out.println(dist.getPDF(i));
		}
		Charter.launchChart(Charter.getChart(range, dist, dist2), "cosine");	
	}
	
	public void testGaussian() {
		// test that we're within 10% of a true Gaussian
		float prior=7f;
		Distribution g=new CosineGaussian(0, 1, prior);
		float[] xs=new float[] {-3, -2, -1, 0, 1, 2, 3};
		float[] cdfs=new float[] {0.0013499672813147567f, 0.022750062887256395f, 0.15865526383236372f, 0.500000000f, 0.8413447361676363f, 0.9772499371127437f, 0.9986500327186852f};
		float[] pdfs=new float[] {0.0044318484119380075f, 0.05399096651318806f, 0.24197072451914337f, 0.3989422804014327f, 0.24197072451914337f, 0.05399096651318806f, 0.0044318484119380075f};
		for (int i=0; i<xs.length; i++) {
			assertEquals(cdfs[i], g.getCDF(xs[i]), 0.1f);
			assertEquals(pdfs[i], g.getPDF(xs[i]), 0.1f);
			assertEquals(pdfs[i]*prior, g.getProbability(xs[i]), prior*0.1f);
		}
	}
}
