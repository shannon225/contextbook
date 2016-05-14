package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import junit.framework.TestCase;

public class LibraryEntryTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(10), new MassTolerance(10), DigestionEnzyme.getEnzyme("trypsin"));
	public void testReverse() {
		double[] massArray = new double[] { 98.06063, 175.11955, 227.10323,
				304.16214, 324.15599, 333, 419.18908, 444, 505.20367,
				532.27314, 555, 618.28773, 650.407259, 666, 713.32082,
				733.31467, 777, 779.449849, 810.37359, 862.35727, 876.502609,
				888, 939.41618, 1018.45838, 1036.46894 };
		float[] intensityArray = new float[] { 1f, 2f, 3f, 4f, 5f, 6f,
				7f, 8f, 9f, 10f, 11f, 12f, 13f, 14f, 15f, 16f, 17f, 18f, 19f,
				20f, 21f, 22f, 23f, 24f, 25f };
		
		LibraryEntry entry=new LibraryEntry("", new HashSet<String>(), 518.73841, (byte)2, "PEPT[+80]IDER", 1, 0.0f, 0.0f, massArray, intensityArray);
		ReverseLibraryEntry reverse=entry.getDecoy(PARAMETERS, false);
		System.out.println(reverse.getPeptideModSeq());;
		assertEquals("EDIT[+80.0]PEPR", reverse.getPeptideModSeq());
		
		double[] reverseMasses=reverse.getMassArray();
		double[] expectedReversedMasses=new double[] { 130.05045900000002, 175.11955, 245.07740900000002, 272.172311, 333.0, 358.161469, 401.21490099999994, 444.0, 498.267661, 538.2850169999999,
				555.0, 636.2619090000001, 650.407259, 666.0, 694.422509, 764.3803670000001, 777.0, 779.449849, 792.399411, 862.3572700000001, 876.502609, 888.0, 938.492048, 1018.4583800000003,
				1036.4689399999997 };
		for (int i = 0; i < reverseMasses.length; i++) {
			assertEquals(expectedReversedMasses[i], reverseMasses[i], 0.001);
		}
	}
}
