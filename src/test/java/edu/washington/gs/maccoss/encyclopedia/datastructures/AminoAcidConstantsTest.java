package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.TestCase;

public class AminoAcidConstantsTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"), false);
	
	public void testGetMass() {
		String pep="AAGPLLTDEC[+57.0]R";
		System.out.println(PARAMETERS.getAAConstants().getChargedMass(pep, (byte)2));
		System.out.println(PARAMETERS.getAAConstants().getMass('C'));
		System.out.println(PARAMETERS.getAAConstants().getMass(pep));
	}
	
	public void testGetChargedMass() {
		AminoAcidConstants NO_MODS = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		String pep="MAREC[+57.0214635]YSPEPTIDESK";
		assertEquals(956.9244, NO_MODS.getChargedMass(pep, (byte)2), 0.0001);
		
		AminoAcidConstants C_MODS = new AminoAcidConstants();
		assertEquals(956.9244, C_MODS.getChargedMass(pep, (byte)2), 0.0001);
	}
}
