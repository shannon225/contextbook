package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Swath;
import junit.framework.TestCase;

public class DotProductTest extends TestCase {
	
	public void testEmptyDotProduct() {
		LibraryEntry entry=getEntry(new double[] {}, new float[] {});
		Swath spectrum=getSwath(new double[] {}, new float[] {});
		assertEquals(0.0f, new DotProduct(new MassTolerance(100.0f)).score(entry, spectrum));
	}
	public void testDotProduct() {
		LibraryEntry entry=getEntry(new double[] {1.0, 29.0, 300.01, 1000.0, 1200.0}, new float[] {7, 7, 2, 3, 7});
		Swath spectrum=getSwath(new double[] {30.0, 300.0, 1001.0, 1300.0}, new float[] {7, 10, 4, 7});
		assertEquals(32.0f, new DotProduct(new MassTolerance(1000.0f)).score(entry, spectrum));
		assertEquals(20.0f, new DotProduct(new MassTolerance(100.0f)).score(entry, spectrum));
		assertEquals(0.0f, new DotProduct(new MassTolerance(10.0f)).score(entry, spectrum));
	}
	
	public LibraryEntry getEntry(double[] masses, float[] intensities) {
		return new LibraryEntry(1, 1, "", 1, 1, 1, masses, intensities);
	}
	public Swath getSwath(double[] masses, float[] intensities) {
		return new Swath("", "", 1, 1, 1, 1, masses, intensities);
	}
}
