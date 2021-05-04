package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassErrorUnitType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import junit.framework.TestCase;

public class MassToleranceTest extends TestCase {
	public static final MassTolerance TOLERANCE=new MassTolerance(10.0);
	public static final MassTolerance AMU_TOLERANCE=new MassTolerance(0.5, MassErrorUnitType.AMU);
	
	public void testNegative() {
		assertTrue(TOLERANCE.equals(-17.0, -17.0));
	}
	
	public void testResolution() {
		MassTolerance resolution=new MassTolerance(30000, MassErrorUnitType.RESOLUTION);
		assertEquals(16.66666666667, resolution.getPpmTolerance(), 0.0000001);
		assertEquals("30000 Resolution", resolution.toString());
	}
	
	public void testMassTolerance() {
		assertEquals(0, TOLERANCE.compareTo(1000.0, 1000.0));
		assertEquals(-1, TOLERANCE.compareTo(1000.0, 1001.0));
		assertEquals(-1, TOLERANCE.compareTo(1000.0, 1000.1));
		assertEquals(0, TOLERANCE.compareTo(1000.0, 1000.01));
		assertEquals(0, TOLERANCE.compareTo(1000.0, 1000.001));
		assertEquals(1, TOLERANCE.compareTo(1000.0, 999.0));
		assertEquals(1, TOLERANCE.compareTo(1000.0, 999.9));
		assertEquals(0, TOLERANCE.compareTo(1000.0, 999.99));
		assertEquals(0, TOLERANCE.compareTo(1000.0, 999.999));

		assertEquals(0, TOLERANCE.compareTo(100.0, 100.0));
		assertEquals(-1, TOLERANCE.compareTo(100.0, 101.0));
		assertEquals(-1, TOLERANCE.compareTo(100.0, 100.1));
		assertEquals(-1, TOLERANCE.compareTo(100.0, 100.01));
		assertEquals(0, TOLERANCE.compareTo(100.0, 100.001));
		assertEquals(1, TOLERANCE.compareTo(100.0, 99.0));
		assertEquals(1, TOLERANCE.compareTo(100.0, 99.9));
		assertEquals(1, TOLERANCE.compareTo(100.0, 99.99));
		assertEquals(0, TOLERANCE.compareTo(100.0, 99.999));

		assertEquals(true, TOLERANCE.equals(1000.0, 1000.0));
		assertEquals(false, TOLERANCE.equals(1000.0, 1001.0));
		assertEquals(false, TOLERANCE.equals(1000.0, 1000.1));
		assertEquals(true, TOLERANCE.equals(1000.0, 1000.01));
		assertEquals(true, TOLERANCE.equals(1000.0, 1000.001));
		assertEquals(false, TOLERANCE.equals(1000.0, 999.0));
		assertEquals(false, TOLERANCE.equals(1000.0, 999.9));
		assertEquals(true, TOLERANCE.equals(1000.0, 999.99));
		assertEquals(true, TOLERANCE.equals(1000.0, 999.999));
	}
	
	public void testGetMatch() {
		double[] values=new double[] {1.0, 1.1, 1.3, 1.399999, 1.4, 1.4000001, 1.5, 2.0};
		assertEquals(new Double(1.4), TOLERANCE.getMatch(values, 1.4).get());
		assertEquals(new Double(1.4000001), TOLERANCE.getMatch(values, 1.40001).get());
	}
	
	public void testGetMatches() {
		double[] values=new double[] {1.0, 1.1, 1.3, 1.399999, 1.4, 1.4000001, 1.5, 2.0};
		double[] matches=TOLERANCE.getMatches(values, 1.4);
		double[] expected=new double[] {1.399999, 1.4, 1.4000001};
		for (int i=0; i<matches.length; i++) {
			assertEquals(expected[i], matches[i]);
		}
	}
	
	public void testAMUMassTolerance() {
		assertEquals(0, AMU_TOLERANCE.compareTo(1000.0, 1000.0));
		assertEquals(-1, AMU_TOLERANCE.compareTo(1000.0, 1001.8));
		assertEquals(1, AMU_TOLERANCE.compareTo(1000.0, 999.2));
		assertEquals(0, AMU_TOLERANCE.compareTo(1000.0, 1000.1));
		assertEquals(0, AMU_TOLERANCE.compareTo(1000.0, 999.9));
		assertEquals(0, AMU_TOLERANCE.compareTo(1000.0, 1000.4));
		assertEquals(0, AMU_TOLERANCE.compareTo(1000.0, 999.6));

		assertEquals(0, AMU_TOLERANCE.compareTo(100.0, 100.0));
		assertEquals(-1, AMU_TOLERANCE.compareTo(100.0, 101.8));
		assertEquals(1, AMU_TOLERANCE.compareTo(100.0, 99.2));
		assertEquals(0, AMU_TOLERANCE.compareTo(100.0, 100.1));
		assertEquals(0, AMU_TOLERANCE.compareTo(100.0, 99.9));
		assertEquals(0, AMU_TOLERANCE.compareTo(100.0, 100.4));
		assertEquals(0, AMU_TOLERANCE.compareTo(100.0, 99.6));
	}
}
