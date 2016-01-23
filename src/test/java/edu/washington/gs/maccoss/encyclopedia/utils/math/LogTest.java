package edu.washington.gs.maccoss.encyclopedia.utils.math;

import junit.framework.TestCase;

public class LogTest extends TestCase {
	public void testLogFactorial() {
		long factorial=1;
		for (int i=1; i<21; i++) {
			factorial=factorial*i;
			float log=Log.log10(factorial);
			
			assertEquals(log, Log.logFactorial(i), 0.0001f);
		}
	}
	
	public void testLog() {
		float dotProduct=1.2046795E+38f;
		int count=5;
		float f=Log.log10(dotProduct)+Log.logFactorial(count);
		assertEquals(40.160053, f, 0.0001);
	}
}
