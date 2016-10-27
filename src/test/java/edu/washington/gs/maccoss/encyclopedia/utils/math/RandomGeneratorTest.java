package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.awt.Color;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import gnu.trove.set.hash.TFloatHashSet;
import junit.framework.TestCase;

public class RandomGeneratorTest extends TestCase {
	public void testRandomSequence() {
		int seed=1;
		String sequence=RandomGenerator.randomSequence(seed);
		assertEquals("CQGLCIDFNMQPEYGIHIIKC", sequence);
	}
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
	
	public void testRandomColor() {
		FragmentIon f1=new FragmentIon(871.0, (byte)4, IonType.y);
		FragmentIon f2=new FragmentIon(971.1, (byte)4, IonType.y);
		Color c1=RandomGenerator.randomColor(f1.toString().hashCode());
		Color c2=RandomGenerator.randomColor(f2.toString().hashCode());
		
		assertEquals(c1.hashCode(), c2.hashCode());
	}
}
