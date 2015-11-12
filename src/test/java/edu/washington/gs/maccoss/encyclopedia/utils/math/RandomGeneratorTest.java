package edu.washington.gs.maccoss.encyclopedia.utils.math;

import gnu.trove.set.hash.TFloatHashSet;
import junit.framework.TestCase;

public class RandomGeneratorTest extends TestCase {
	public void testRandomInt() {
		TFloatHashSet set=new TFloatHashSet();
		int[] count=new int[1000];
		for (int i=0; i<100000; i++) {
			float r=RandomGenerator.random(i);
			// asserts you never get the same number in seeds of 0 to 100,000
			assertFalse(set.contains(r));
			set.add(r);
			int index=(int)(count.length*r);
			count[index]++;
		}
		
		// asserts regularity
		for (int i=0; i<count.length; i++) {
			assertTrue(count[i]>95);
			assertTrue(count[i]<105);
		}
	}
}
