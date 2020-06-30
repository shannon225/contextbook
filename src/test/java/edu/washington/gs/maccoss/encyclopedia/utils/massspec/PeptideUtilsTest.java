package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.TestCase;

public class PeptideUtilsTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"), false, true, false);
	
	public void testGetExpectedChargeState() {
		assertEquals(1, PeptideUtils.getExpectedChargeState("LACDEFQFEDCAI"));
		assertEquals(2, PeptideUtils.getExpectedChargeState("LACDEFQFEDCAIR"));
		assertEquals(3, PeptideUtils.getExpectedChargeState("LACDEFQKEDCAIR"));
		assertEquals(3, PeptideUtils.getExpectedChargeState("LACHEFQFEDCAIR"));
		assertEquals(4, PeptideUtils.getExpectedChargeState("LACHEFQFEDRAIR"));
		assertEquals(5, PeptideUtils.getExpectedChargeState("LRCHEFQFEDRAIR"));
	}
	
	public void testGetModIndicies() {
		String sequence="A[+42.0]QRHS[+79.96633]DSCCSLEEK";
		assertEquals(1, PeptideUtils.getNumberOfMods(sequence, 80));
		assertEquals(5, PeptideUtils.getModIndicies(sequence, 80)[0]);
	}
	
	public void testDecoys() {
		HashSet<String> backgroundProteome=new HashSet<String>();
		String seq="LLACDEFQFEDCAIR";
		backgroundProteome.add(seq);
		
		// dumb decoy method sees these as different peptides
		assertEquals("LIACDEFQFEDCALR", PeptideUtils.getDecoy(seq, backgroundProteome, PARAMETERS));
		// but they share too many ions, so actually must shuffle
		assertEquals("LCAEFCILFAQDEDR", PeptideUtils.getSmartDecoy(seq, (byte)2, backgroundProteome, PARAMETERS));
		

		assertEquals("LCAEFCILFAQDEDR", PeptideUtils.getSmartDecoy(seq, (byte)3, backgroundProteome, PARAMETERS));
	}

	public void testReverse() {
		String s=PeptideUtils.reverse("ABC[+57]DEFGHIJK", PARAMETERS.getAAConstants());
		assertEquals("AJIHGFEDC[+57.0214635]BK", s);
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
		String sequence="A[+42.0]QRHS[+79.96633]DSC[+57.0214635]C[+57.0214635]SLEEK";
		String peptideModSeq=PeptideUtils.formatForSkyline(sequence);
		assertEquals("A[+42.0]QRHS[+80.0]DSC[+57.0]C[+57.0]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForSkylinePeakBoundaries(sequence);
		assertEquals("A[+42]QRHS[+80]DSC[+57]C[+57]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForEncyclopeDIA(sequence);
		assertEquals("A[+42.0]QRHS[+79.96633]DSC[+57.0214635]C[+57.0214635]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForEncyclopeDIA("A[+42.0]QRHS[+79.96633]DSC[+57.0]C[+57.0]SLEEK");
		assertEquals("A[+42.0]QRHS[+79.96633]DSC[+57.0]C[+57.0]SLEEK", peptideModSeq);

		sequence="Q[-17.0]QRHS[+79.96633]DSC[+57.0214635]C[+57.0214635]SLEEK";
		peptideModSeq=PeptideUtils.formatForSkyline(sequence);
		assertEquals("Q[-17.0]QRHS[+80.0]DSC[+57.0]C[+57.0]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForSkylinePeakBoundaries(sequence);
		assertEquals("Q[-17]QRHS[+80]DSC[+57]C[+57]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForEncyclopeDIA(sequence);
		assertEquals("Q[-17.0]QRHS[+79.96633]DSC[+57.0214635]C[+57.0214635]SLEEK", peptideModSeq);

		peptideModSeq=PeptideUtils.formatForEncyclopeDIA("Q[-17.0]QRHS[+79.96633]DSC[+57.0]C[+57.0]SLEEK");
		assertEquals("Q[-17.0]QRHS[+79.96633]DSC[+57.0]C[+57.0]SLEEK", peptideModSeq);
	}
	
	public void testGetCorrectedMasses() {
		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		String sequence="A[+42.0]QRHS[+79.96633]DSCCSLEEK";
		assertEquals("A[+42.010565]QRHS[+79.966331]DSCCSLEEK", PeptideUtils.getCorrectedMasses(sequence, aaConstants));
		sequence="A[+42.0]QRHS[+80.0]DSCCSLEEK";
		assertEquals("A[+42.010565]QRHS[+79.966331]DSCCSLEEK", PeptideUtils.getCorrectedMasses(sequence, aaConstants));
		sequence="Q[-17]QRHS[+80.0]DSCCSLEEK";
		assertEquals("Q[-17.026549]QRHS[+79.966331]DSCCSLEEK", PeptideUtils.getCorrectedMasses(sequence, aaConstants));
		sequence="Q[-17.026549]QRHS[+79.96633]DSCCSLEEK";
		assertEquals("Q[-17.026549]QRHS[+79.966331]DSCCSLEEK", PeptideUtils.getCorrectedMasses(sequence, aaConstants));
	}
	
	public void testGetAAs() {
		AminoAcidConstants aminoAcidConstants = AminoAcidConstants.createEmptyFixedAndVariable();
		for (int i=0; i<100; i++) {
			StringBuilder sb=new StringBuilder();
			double length=Math.random()*30+5;
			for (int j=0; j<length; j++) {
				int index=(int)(AminoAcidConstants.AAs.length*Math.random());
				if (index>=AminoAcidConstants.AAs.length) index=AminoAcidConstants.AAs.length-1;
				sb.append(AminoAcidConstants.AAs[index]);
				if (Math.random()>0.8) {
					sb.append("[+"+(100+Math.random()*100)+"]"); // outside any mass reformatting
				}
			}
			String sequence=sb.toString();
			String processedSequence=General.toString(PeptideUtils.getAAs(sequence, aminoAcidConstants), "");
			assertEquals(sequence, processedSequence);
		}
	}
	
	public void testGetMasses() {
		String sequence="PEPTIDER";
		double[] expected=new double[] {97.0528, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011};
		double[] masses=PeptideUtils.getPeptideModel(sequence, PARAMETERS.getAAConstants()).getMasses();
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.1);
		}

		sequence="PEPT[+80]IDER";
		expected=new double[] {97.0528, 129.0426, 97.0528, 101.0477+80.0, 113.0841, 115.027, 129.0426, 156.1011};
		masses=PeptideUtils.getPeptideModel(sequence, PARAMETERS.getAAConstants()).getMasses();
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.1);
		}

		sequence="PE[-17]PTIDER";
		expected=new double[] {97.0528, 129.0426-17.0, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011};
		masses=PeptideUtils.getPeptideModel(sequence, PARAMETERS.getAAConstants()).getMasses();
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.1);
		}

		sequence="[-17]PEPTIDER";
		expected=new double[] {97.0528-17.0, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011};
		masses=PeptideUtils.getPeptideModel(sequence, PARAMETERS.getAAConstants()).getMasses();
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.1);
		}

		sequence="[+42]PEPTIDER";
		expected=new double[] {97.0528+42.0, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011};
		masses=PeptideUtils.getPeptideModel(sequence, PARAMETERS.getAAConstants()).getMasses();
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.1);
		}

		sequence="PEPTIDER[+14]";
		expected=new double[] {97.0528, 129.0426, 97.0528, 101.0477, 113.0841, 115.027, 129.0426, 156.1011+14.0};
		masses=PeptideUtils.getPeptideModel(sequence, PARAMETERS.getAAConstants()).getMasses();
		for (int i=0; i<masses.length; i++) {
			assertEquals(expected[i], masses[i], 0.1);
		}
	}
}
