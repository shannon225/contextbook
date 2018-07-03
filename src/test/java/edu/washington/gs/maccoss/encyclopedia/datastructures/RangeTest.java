package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;

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

	public void testWidestRange() {
		ArrayList<Range> ranges=new ArrayList<Range>();
		ranges.add(new Range(10f, 20f));
		ranges.add(new Range(30f, 40f));
		ranges.add(new Range(-10f, 20f));
		ranges.add(new Range(10f, 20f));
		
		Range widest=Range.getWidestRange(ranges);
		assertEquals(-10f, widest.getStart());
		assertEquals(40f, widest.getStop());
	}
	
	public void testContains() {
		Range r1=new Range(100, 200);
		Range r2=new Range(199, 199);
		Range r3=new Range(199, 201);
		Range r4=new Range(200, 300);
		Range r5=new Range(300, 400);

		Range r6=new Range(101, 101);
		Range r7=new Range(99, 101);
		Range r8=new Range(98, 99);
		
		Range r9=new Range(0, 300);
		
		assertTrue(r1.contains(r2));
		assertTrue(r1.contains(r3));
		assertTrue(r1.contains(r4));
		assertFalse(r1.contains(r5));
		assertTrue(r1.contains(r6));
		assertTrue(r1.contains(r7));
		assertFalse(r1.contains(r8));
		assertTrue(r1.contains(r9));
	}
}
