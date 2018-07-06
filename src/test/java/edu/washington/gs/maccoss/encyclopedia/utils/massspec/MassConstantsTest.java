package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import junit.framework.TestCase;

public class MassConstantsTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
	
	public void testGetMass() {
		String sequence="PEPTIDER";
		assertEquals(955.46112, PARAMETERS.getAAConstants().getMass(sequence)+18.01042, 0.001);
	}
	
	public void testGetNeutralLoss() {
		AminoAcidConstants constants = new AminoAcidConstants();
		assertTrue(constants.getNeutralLoss('S', 80.0)>0.0);
		assertTrue(constants.getNeutralLoss('S', 80.1)==0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331)>0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331+0.00001)>0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331+0.0001)==0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331-0.00001)>0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331-0.0001)==0.0);
	}
}
