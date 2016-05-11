package edu.washington.gs.maccoss.encyclopedia.datastructures;

import junit.framework.TestCase;

public class RangeTest extends TestCase {
	public void testInterpolation() {
		Range r=new Range(30f, 40f);
		// middle
		assertEquals(500f, r.linearInterp(35f, 0f, 1000f));
		assertEquals(45f, r.linearInterp(35f, 40f, 50f));
		
		// boundaries
		assertEquals(40f, r.linearInterp(30f, 40f, 50f));
		assertEquals(50f, r.linearInterp(40f, 40f, 50f));
		

		// middle
		assertEquals(35f, r.mapBackToRange(500f, 0f, 1000f));
		assertEquals(35f, r.mapBackToRange(45f, 40f, 50f));
		
		// boundaries
		assertEquals(30f, r.mapBackToRange(40f, 40f, 50f));
		assertEquals(40f, r.mapBackToRange(50f, 40f, 50f));
	}

}
