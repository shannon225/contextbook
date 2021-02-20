package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.filereaders.PTMMap.PostTranslationalModification;
import junit.framework.TestCase;

public class PTMMapTest extends TestCase {

	public void testGetPTM() {
		assertEquals("Carbamidomethyl", PTMMap.getPTM("Carbamidomethyl", "C").getName());
		assertEquals("Carbamidomethyl", PTMMap.getPTM("UNIMOD:4").getName());
		assertEquals(PostTranslationalModification.nothing.getName(), PTMMap.getPTM("UniMod:4").getName());
		assertEquals("Carbamidomethyl", PTMMap.getPTM(57.0214635, "C").getName());
		assertEquals("Carbamidomethyl", PTMMap.getPTM(57.021463, "C").getName());
		assertEquals("Carbamidomethyl", PTMMap.getPTM(57.02146, "C").getName());
		
		assertFalse("Carbamidomethyl".equals(PTMMap.getPTM(57.0214, "C").getName())); // failure at float
		assertFalse("Carbamidomethyl".equals(PTMMap.getPTM(15.994915, "C").getName())); // wrong mass

		assertEquals("Oxidation", (PTMMap.getPTM(15.994915, "MW").getName()));
		assertEquals("Oxidation", (PTMMap.getPTM(15.994915, "M").getName()));
		assertEquals("Oxidation", (PTMMap.getPTM(15.994915, "W").getName()));
		assertFalse("Oxidation".equals(PTMMap.getPTM(15.994915, "I").getName())); // wrong AA
	}
}
