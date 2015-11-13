package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import junit.framework.TestCase;

public class MassConstantsTest extends TestCase {
	public void testGetMass() {
		String sequence="PEPTIDER";
		assertEquals(955.46112, MassConstants.getMass(sequence)+18.01042, 0.001);
	}
	
	public void testGetNeutralLoss() {
		assertTrue(MassConstants.getNeutralLoss(80.0)>0.0);
		assertTrue(MassConstants.getNeutralLoss(80.1)==0.0);
		assertTrue(MassConstants.getNeutralLoss(79.966331)>0.0);
		assertTrue(MassConstants.getNeutralLoss(79.966331+0.00001)>0.0);
		assertTrue(MassConstants.getNeutralLoss(79.966331+0.0001)==0.0);
		assertTrue(MassConstants.getNeutralLoss(79.966331-0.00001)>0.0);
		assertTrue(MassConstants.getNeutralLoss(79.966331-0.0001)==0.0);
	}
}
