package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import junit.framework.TestCase;

public class DotProductTest extends TestCase {
	
	public void testEmptyDotProduct() {
		LibraryEntry entry=getEntry(new double[] {}, new float[] {});
		Stripe spectrum=getStripe(new double[] {}, new float[] {});
		assertEquals(0.0f, new DotProduct(new MassTolerance(100.0f)).score(entry, spectrum, new float[] {}, null));
	}
	public void testDotProduct() {
		LibraryEntry entry=getEntry(new double[] {1.0, 29.0, 300.01, 1000.0, 1200.0}, new float[] {7, 7, 2, 3, 7});
		Stripe spectrum=getStripe(new double[] {30.0, 300.0, 1001.0, 1300.0}, new float[] {7, 10, 4, 7});
		assertEquals(32.0f, new DotProduct(new MassTolerance(1000.0f)).score(entry, spectrum, new float[] {}, null));
		assertEquals(20.0f, new DotProduct(new MassTolerance(100.0f)).score(entry, spectrum, new float[] {}, null));
		assertEquals(0.0f, new DotProduct(new MassTolerance(10.0f)).score(entry, spectrum, new float[] {}, null));
	}
	
	public LibraryEntry getEntry(double[] masses, float[] intensities) {
		return new LibraryEntry(1, (byte)1, "", 1, 1, 1, masses, intensities);
	}
	public Stripe getStripe(double[] masses, float[] intensities) {
		return new Stripe("", "", 1, 1, 1, 1, masses, intensities);
	}
}
