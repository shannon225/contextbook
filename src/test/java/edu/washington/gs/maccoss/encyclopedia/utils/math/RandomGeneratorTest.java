package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.awt.Color;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import gnu.trove.set.hash.TFloatHashSet;
import junit.framework.TestCase;

public class RandomGeneratorTest extends TestCase {
	public void testExpectedRandomNumbers() {
		assertEquals(12345, RandomGenerator.randomInt(0));
		assertEquals(1103527590, RandomGenerator.randomInt(1));
		assertEquals(-1103502900, RandomGenerator.randomInt(-1));
		assertEquals(1043980748, RandomGenerator.randomInt(Integer.MAX_VALUE));
		assertEquals(-2147471303, RandomGenerator.randomInt(Integer.MIN_VALUE));
	}
	
	public void testRandomSequence() {
		int seed=1;
		String sequence=RandomGenerator.randomSequence(seed);
		assertEquals("CQGLCIDFNMQPEYGIHIIKC", sequence);
	}
	
	/**
	 * not quite as regular as the ANSI C method, but still pretty good
	 */
	public void testRandomIntAlt() {
		TFloatHashSet set=new TFloatHashSet();
		int[] count=new int[1000];
		int seed=1;
		for (int i=0; i<100000; i++) {
			seed=RandomGenerator.randomInt(seed);
			float r=RandomGenerator.floatFromRandomInt(RandomGenerator.randomIntAlt(i));
			// asserts you never get the same number in seeds of 0 to 100,000
			assertFalse(set.contains(r));
			set.add(r);
			int index=(int)(count.length*r);
			count[index]++;
		}
		
		// asserts regularity
		for (int i=0; i<count.length; i++) {
			assertTrue(count[i]>85);
			assertTrue(count[i]<115);
		}
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
		Ion f1=new FragmentIon(871.0, (byte)4, IonType.y);
		Ion f2=new FragmentIon(971.1, (byte)4, IonType.y);
		Color c1=RandomGenerator.randomColor(f1.toString().hashCode());
		Color c2=RandomGenerator.randomColor(f2.toString().hashCode());
		
		assertEquals(c1.hashCode(), c2.hashCode());
	}
}
