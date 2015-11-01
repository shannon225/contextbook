package edu.washington.gs.maccoss.encyclopedia.algorithms;

import junit.framework.TestCase;

public class MassToleranceTest extends TestCase {
	public void testMassTolerance() {
		MassTolerance t=new MassTolerance(10.0);
		assertEquals(0, t.compareTo(1000.0, 1000.0));
		assertEquals(-1, t.compareTo(1000.0, 1001.0));
		assertEquals(-1, t.compareTo(1000.0, 1000.1));
		assertEquals(0, t.compareTo(1000.0, 1000.01));
		assertEquals(0, t.compareTo(1000.0, 1000.001));
		assertEquals(1, t.compareTo(1000.0, 999.0));
		assertEquals(1, t.compareTo(1000.0, 999.9));
		assertEquals(0, t.compareTo(1000.0, 999.99));
		assertEquals(0, t.compareTo(1000.0, 999.999));

		assertEquals(true, t.equals(1000.0, 1000.0));
		assertEquals(false, t.equals(1000.0, 1001.0));
		assertEquals(false, t.equals(1000.0, 1000.1));
		assertEquals(true, t.equals(1000.0, 1000.01));
		assertEquals(true, t.equals(1000.0, 1000.001));
		assertEquals(false, t.equals(1000.0, 999.0));
		assertEquals(false, t.equals(1000.0, 999.9));
		assertEquals(true, t.equals(1000.0, 999.99));
		assertEquals(true, t.equals(1000.0, 999.999));
	}
}
