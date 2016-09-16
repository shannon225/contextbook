package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import junit.framework.TestCase;

public class PeptideUtilsTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"));
	
	public void testDecoys() {
		HashSet<String> backgroundProteome=new HashSet<String>();
		String seq="LACDEFQFEDCAIR";
		backgroundProteome.add(seq);
		
		// dumb decoy method sees these as different peptides
		assertEquals("IACDEFQFEDCALR", PeptideUtils.getDecoy(seq, backgroundProteome, PARAMETERS));
		// but they share too many ions, so actually must shuffle
		assertEquals("QFIDDACECEFALR", PeptideUtils.getSmartDecoy(seq, (byte)2, backgroundProteome, PARAMETERS));
		

		assertEquals("QFIDDACECEFALR", PeptideUtils.getSmartDecoy(seq, (byte)3, backgroundProteome, PARAMETERS));
	}

	public void testReverse() {
		String s=PeptideUtils.reverse("ABC[+57]DEFGHIJK", PARAMETERS);
		assertEquals("JIHGFEDC[+57.0]BAK", s);
	}
	
	public void testShuffle() {
		String s="ABC[+57]DEFGHIJK";
		 
		HashSet<String> set=new HashSet<String>();
		set.add(s);
		for (int i=0; i<1000; i++) {
			String shuffle=PeptideUtils.shuffle(s, PARAMETERS);
			s=PeptideUtils.shuffle(s, PARAMETERS);
			// asserts random from the same seed always returns the same sequence
			assertEquals(shuffle, s);

			// asserts random 11mer is never reused in 1000 sequences
			assertFalse(set.contains(s));
			set.add(s);
		}
	}
	
	public void testSkylinePeptideModSeq() {
		String sequence="A[+42.0]QRHS[+79.96633]DSCCSLEEK";
		String peptideModSeq=PeptideUtils.formatForSkyline(sequence, PARAMETERS.getAAConstants());
		assertEquals("A[+42.0]QRHS[+80.0]DSC[+57.0]C[+57.0]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForSkylinePeakBoundaries(sequence, PARAMETERS.getAAConstants());
		assertEquals("A[+42]QRHS[+80]DSC[+57]C[+57]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForEncyclopeDIA(sequence, PARAMETERS.getAAConstants());
		assertEquals("A[+42.0]QRHS[+79.96633]DSC[+57.0214635]C[+57.0214635]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForEncyclopeDIA("A[+42.0]QRHS[+79.96633]DSC[+57.0]C[+57.0]SLEEK", PARAMETERS.getAAConstants());
		assertEquals("A[+42.0]QRHS[+79.96633]DSC[+57.0]C[+57.0]SLEEK", peptideModSeq);

		sequence="Q[-17.0]QRHS[+79.96633]DSCCSLEEK";
		peptideModSeq=PeptideUtils.formatForSkyline(sequence, PARAMETERS.getAAConstants());
		assertEquals("Q[-17.0]QRHS[+80.0]DSC[+57.0]C[+57.0]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForSkylinePeakBoundaries(sequence, PARAMETERS.getAAConstants());
		assertEquals("Q[-17]QRHS[+80]DSC[+57]C[+57]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForEncyclopeDIA(sequence, PARAMETERS.getAAConstants());
		assertEquals("Q[-17.0]QRHS[+79.96633]DSC[+57.0214635]C[+57.0214635]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForEncyclopeDIA("Q[-17.0]QRHS[+79.96633]DSC[+57.0]C[+57.0]SLEEK", PARAMETERS.getAAConstants());
		assertEquals("Q[-17.0]QRHS[+79.96633]DSC[+57.0]C[+57.0]SLEEK", peptideModSeq);
	}
	
	public void testGetMasses() {
		String sequence="PEPTIDER";
		double[] expected=new double[] {97.0528, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011};
		double[] masses=PeptideUtils.getMasses(sequence, PARAMETERS.getAAConstants()).x;
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.001);
		}

		sequence="PEPT[+80]IDER";
		expected=new double[] {97.0528, 129.0426, 97.0528, 101.0477+80.0, 113.0841, 115.027, 129.0426, 156.1011};
		masses=PeptideUtils.getMasses(sequence, PARAMETERS.getAAConstants()).x;
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.001);
		}

		sequence="PE[-17]PTIDER";
		expected=new double[] {97.0528, 129.0426-17.0, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011};
		masses=PeptideUtils.getMasses(sequence, PARAMETERS.getAAConstants()).x;
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.001);
		}

		sequence="[-17]PEPTIDER";
		expected=new double[] {97.0528-17.0, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011};
		masses=PeptideUtils.getMasses(sequence, PARAMETERS.getAAConstants()).x;
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.001);
		}

		sequence="[+42]PEPTIDER";
		expected=new double[] {97.0528+42.0, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011};
		masses=PeptideUtils.getMasses(sequence, PARAMETERS.getAAConstants()).x;
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.001);
		}

		sequence="PEPTIDER[+14]";
		expected=new double[] {97.0528, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011+14.0};
		masses=PeptideUtils.getMasses(sequence, PARAMETERS.getAAConstants()).x;
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.001);
		}
	}
}
