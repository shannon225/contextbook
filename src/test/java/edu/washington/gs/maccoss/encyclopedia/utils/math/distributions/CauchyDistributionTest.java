package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import junit.framework.TestCase;

public class CauchyDistributionTest extends TestCase {
	public static void main(String[] args) {
		CauchyDistribution dist=new CauchyDistribution(0, 0.2f, 5);
		Range range=new Range(-5, 5);
		Charter.launchChart(Charter.getChart(range, dist), "CauchyDistribution");
	}

	public void testGaussian() {
		// test that we're within 10% of a true Gaussian
		CauchyDistribution g=new CauchyDistribution(0, 0.2f, 5);
		float[] xs=new float[] {-3, -2, -1, 0, 1, 2, 3};
		float[] cdfs=new float[] {0.02118930494975968f, 0.03172551790017597f, 0.0628329591011525f, 0.5f, 0.9371670408988475f, 0.968274482099824f, 0.9788106950502403f};
		float[] pdfs=new float[] {0.007042254223119523f, 0.015757915387774907f, 0.061213440492714814f, 1.591549407203019f, 0.061213440492714814f, 0.015757915387774907f, 0.007042254223119523f};
		for (int i=0; i<xs.length; i++) {
			assertEquals(cdfs[i], g.getCDF(xs[i]), 0.0001f);
			assertEquals(pdfs[i], g.getPDF(xs[i]), 0.0001f);
		}
	}
}
