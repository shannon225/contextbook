package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.filereaders.PTMMap.PostTranslationalModification;
import junit.framework.TestCase;

public class PTMMapTest extends TestCase {

	public void testGetPTM() {
		assertEquals("Carbamidomethyl", PTMMap.getPTM("Carbamidomethyl", "C").getName());
		assertEquals("Carbamidomethyl", PTMMap.getPTM("UNIMOD:4").getName());
		assertEquals(PostTranslationalModification.nothing.getName(), PTMMap.getPTM("UniMod:4").getName());
	}
}
