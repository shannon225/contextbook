package edu.washington.gs.maccoss.encyclopedia.utils.math;

import junit.framework.TestCase;

public class GeneralTest extends TestCase {
	public void testConcatenate() {
		float[] a1=new float[] {1, 2, 3};
		float[] a2=new float[] {4, 5, 6};
		float[] a3=new float[] {7, 8, 9};

		float[] r=General.concatenate(a1, a2, a3);
		float[] expected=new float[] {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
		for (int i=0; i<r.length; i++) {
			assertEquals(expected[i], r[i]);
		}
	}

}
