package edu.washington.gs.maccoss.encyclopedia.algorithms;

import static org.junit.Assert.*;

import org.junit.Test;

public class SSRCalcTest {

	@Test
	public void testSAVDLSCSR() {
		assertEquals(12.550075369248415, SSRCalc.getHydrophobicity("SAVDLSCSR"), 0.000001);
	}

}
