package edu.washington.gs.maccoss.encyclopedia.datastructures;

import junit.framework.TestCase;

public class ModificationMassMapTest extends TestCase {
	public void testMap() {
		String value="C=57.0";
		ModificationMassMap map=new ModificationMassMap(value);
		assertEquals(57.0, map.getVariableMod('C'));
		assertEquals(ModificationMassMap.MISSING, map.getVariableMod('c'));
		assertEquals(value, map.toString());

		map=new ModificationMassMap("K=8.014199,W=15.994915,nQ=-17.026549,nC=-17.026549,a=42.010565");
		assertEquals(8.014199, map.getVariableMod('K'));
		assertEquals(15.994915, map.getVariableMod('W'));

		assertEquals(ModificationMassMap.MISSING, map.getVariableMod('n'));
		assertEquals(-17.026549, map.getNTermMod('Q'));
		assertEquals(-17.026549, map.getNTermMod('C'));
		assertEquals(ModificationMassMap.MISSING, map.getNTermMod('n'));
		assertEquals(42.010565, map.getProteinNTermMod('M'));
		assertEquals(42.010565, map.getProteinNTermMod('Q'));
		assertEquals("K=8.014199,W=15.994915,a=42.010565,nC=-17.026549,nQ=-17.026549", map.toString());
		
	}

}
