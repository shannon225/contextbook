package edu.washington.gs.maccoss.encyclopedia.datastructures;

import junit.framework.TestCase;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;

public class LibraryEntryTest extends TestCase {
	public void testReverse() {
		double[] massArray = new double[] { 98.06063, 175.11955, 227.10323,
				304.16214, 324.15599, 333, 419.18908, 444, 505.20367,
				532.27314, 555, 618.28773, 650.407259, 666, 713.32082,
				733.31467, 777, 779.449849, 810.37359, 862.35727, 876.502609,
				888, 939.41618, 1018.45838, 1036.46894 };
		float[] intensityArray = new float[] { 48f, 41f, 97f, 22f, 5f, 20f,
				68f, 74f, 70f, 87f, 4f, 50f, 51f, 90f, 73f, 66f, 12f, 9f, 76f,
				15f, 34f, 62f, 81f, 94f, 50f };
		
		LibraryEntry entry=new LibraryEntry(518.73841, (byte)2, "PEPT[+80]IDER", 1, 0.0f, 0.0f, massArray, intensityArray);
		ReverseLibraryEntry reverse=entry.getReverse(new MassTolerance(10.0f));
		assertEquals("REDIT[+80]PEP", reverse.getPeptideModSeq());
		
		double[] reverseMasses=reverse.getMassArray();
		double[] expectedReversedMasses = new double[] { 116.0712, 157.10898,
				245.1138, 286.15158, 333, 342.1666, 401.17858, 444, 514.26268,
				523.2143, 555, 636.2984, 650.407259, 666, 695.31038, 751.3254,
				777, 779.449849, 792.36318, 876.502609, 880.368, 888,
				921.40578, 1018.45858, 1036.4691 };
		for (int i = 0; i < reverseMasses.length; i++) {
			assertEquals(expectedReversedMasses[i], reverseMasses[i], 0.001);
		}
	}
}
