package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import junit.framework.TestCase;

public class AminoAcidConstantsTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"));
	
	public void testGetMass() {
		String pep="AAGPLLTDEC[+57.0]R";
		System.out.println(PARAMETERS.getAAConstants().getChargedMass(pep, (byte)2));
		System.out.println(PARAMETERS.getAAConstants().getMass('C'));
		System.out.println(PARAMETERS.getAAConstants().getMass(pep));
	}
}
