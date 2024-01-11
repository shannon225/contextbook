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
		
		sequence="MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
		assertEquals(8560.6238, PARAMETERS.getAAConstants().getMass(sequence)+MassConstants.oh2+MassConstants.protonMass, 0.001);
	}
	
	public void testGetModificationMasses() {
		String[] sequences=new String[] {"EEEAIQLDGLN[+203.079373]ASQIR",
		"ELISN[+203.079373]ASDALDK",
		"HNN[+203.079373]DTQHIWESDSNEFSVIADPR",
		"GVVDSDDLPLN[+203.079373]VSR",
		"YN[+203.079373]DTFWK",
		"LGVIEDHSN[+203.079373]R",
		"LGVIEDHSN[+203.079373]RTR"};
		
		double[] targets=new double[] {1987.9696849999998,
				1477.714638,
				2913.280012,
				1687.8263140000001,
				1175.513352,
				1341.6523129999998,
				1598.8011029999998};
		
		for (int i = 0; i < sequences.length; i++) {
			assertEquals(targets[i], PARAMETERS.getAAConstants().getMass(sequences[i])+18.01042, 0.001);
		}
	}
	
	public void testGetNeutralLoss() {
		AminoAcidConstants constants = new AminoAcidConstants();
		assertTrue(constants.getNeutralLoss('S', 80.0)>0.0);
		assertTrue(constants.getNeutralLoss('S', 80.1)==0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331)>0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331+0.00001)>0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331+0.001)==0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331-0.00001)>0.0);
		assertTrue(constants.getNeutralLoss('S', 79.966331-0.001)==0.0);
	}
}
