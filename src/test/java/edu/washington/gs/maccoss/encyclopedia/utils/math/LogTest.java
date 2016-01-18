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
}
