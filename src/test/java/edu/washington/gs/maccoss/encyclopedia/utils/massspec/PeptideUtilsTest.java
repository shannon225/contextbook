package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
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
		assertEquals("FDFCDQECELAAIR", PeptideUtils.getSmartDecoy(seq, (byte)2, backgroundProteome, PARAMETERS));
		

		assertEquals("FDFCDQECELAAIR", PeptideUtils.getSmartDecoy(seq, (byte)3, backgroundProteome, PARAMETERS));
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
}
