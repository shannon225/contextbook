package edu.washington.gs.maccoss.encyclopedia.utils.io;

import junit.framework.TestCase;

public class VersionTest extends TestCase {
	public void testVersion() {
		Version v1=new Version("1.2");
		Version v2=new Version("1.2.0");
		Version v3=new Version("1.2.1");
		Version v4=new Version("1.1.9");
		Version v5=new Version("0.9.8");
		
		assertEquals("1.2.0", v1.toString());
		assertEquals("1.2.0", v2.toString());
		assertEquals("1.2.1", v3.toString());
		assertEquals("1.1.9", v4.toString());
		assertEquals("0.9.8", v5.toString());
		
		assertTrue(v4.compareTo(v5)==1);
		assertTrue(v5.compareTo(v3)==-1);
		assertTrue(v1.compareTo(v2)==0);
		assertTrue(v1.equals(v2));
		assertFalse(v1.equals(v3));
		
		assertTrue(v3.amIAbove(v2));
		assertFalse(v4.amIAbove(v2));
		assertFalse(v2.amIAbove(v2));
	}

}
