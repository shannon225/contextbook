package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.HashSet;

import junit.framework.TestCase;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;

public class PeptideUtilsTest extends TestCase {
	private static final DigestionEnzyme ENZYME=DigestionEnzyme.getEnzyme("trypsin");
	private static final SearchParameters PARAMETERS=new SearchParameters(FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), ENZYME);
	
	public void testDecoys() {
		HashSet<String> backgroundProteome=new HashSet<String>();
		String seq="LACDEFQFEDCAIR";
		backgroundProteome.add(seq);
		
		// dumb decoy method sees these as different peptides
		assertEquals("IACDEFQFEDCALR", PeptideUtils.getDecoy(seq, ENZYME, backgroundProteome));
		// but they share too many ions, so actually must shuffle
		assertEquals("FDFCDQECELAAIR", PeptideUtils.getSmartDecoy(seq, backgroundProteome, PARAMETERS));
	}

	public void testReverse() {
		String s=PeptideUtils.reverse("ABC[+57]DEFGHIJK", ENZYME);
		assertEquals("JIHGFEDC[+57]BAK", s);
	}
	
	public void testShuffle() {
		String s="ABC[+57]DEFGHIJK";
		 
		HashSet<String> set=new HashSet<String>();
		set.add(s);
		for (int i=0; i<1000; i++) {
			String shuffle=PeptideUtils.shuffle(s, ENZYME);
			s=PeptideUtils.shuffle(s, ENZYME);
			// asserts random from the same seed always returns the same sequence
			assertEquals(shuffle, s);

			// asserts random 11mer is never reused in 1000 sequences
			assertFalse(set.contains(s));
			set.add(s);
		}
		
		
	}
}
