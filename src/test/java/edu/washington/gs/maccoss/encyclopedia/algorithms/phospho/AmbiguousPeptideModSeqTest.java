package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

public class AmbiguousPeptideModSeqTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"));
	
	public void testAmbiguity() {;
		assertEquals("(S[+79.96633])SSR", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.96633]SSR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]SS)R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.96633]SSR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SS[+79.96633])SR", AmbiguousPeptideModSeq.getLeftAmbiguity("SS[+79.96633]SR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.96633]S)R", AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.96633]SR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SSS[+79.96633])R", AmbiguousPeptideModSeq.getLeftAmbiguity("SSS[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("SS(S[+79.96633])R", AmbiguousPeptideModSeq.getRightAmbiguity("SSS[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		

		assertEquals("(S[+79.96633])(S[+79.96633])SR", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.96633]S[+79.96633]SR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633])(S[+79.96633]S)R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.96633]S[+79.96633]SR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SS[+79.96633])(S[+79.96633])R", AmbiguousPeptideModSeq.getLeftAmbiguity("SS[+79.96633]S[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.96633])(S[+79.96633])R", AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.96633]S[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633])(SS[+79.96633])R", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.96633]SS[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]S)(S[+79.96633])R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.96633]SS[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());

		assertEquals("(S[+79.96633])SSR", AmbiguousPeptideModSeq.getUnambigous("S[+79.96633]SSR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.96633])SR", AmbiguousPeptideModSeq.getUnambigous("SS[+79.96633]SR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("SS(S[+79.96633])R", AmbiguousPeptideModSeq.getUnambigous("SSS[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633])(S[+79.96633])SR", AmbiguousPeptideModSeq.getUnambigous("S[+79.96633]S[+79.96633]SR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.96633])(S[+79.96633])R", AmbiguousPeptideModSeq.getUnambigous("SS[+79.96633]S[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633])S(S[+79.96633])R", AmbiguousPeptideModSeq.getUnambigous("S[+79.96633]SS[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());

		assertEquals("(S[+79.96633]SS)R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.96633]SSR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SS[+79.96633]S)R", AmbiguousPeptideModSeq.getFullyAmbiguous("SS[+79.96633]SR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SSS[+79.96633])R", AmbiguousPeptideModSeq.getFullyAmbiguous("SSS[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]S[+79.96633]S)R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.96633]S[+79.96633]SR",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SS[+79.96633]S[+79.96633])R", AmbiguousPeptideModSeq.getFullyAmbiguous("SS[+79.96633]S[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]SS[+79.96633])R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.96633]SS[+79.96633]R",PARAMETERS.getAAConstants()).getPeptideAnnotation());
	}
	
	public void testSets() {
		AmbiguousPeptideModSeq s=AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.96633]SR",PARAMETERS.getAAConstants());
		TIntHashSet[] sets=s.getAmbiguityGroups();
		assertTrue(sets.length==1);
	}
	
	public void testIsLocalized() {
		String s="(S[+80.0])SSSR";
		assertTrue(AmbiguousPeptideModSeq.isLocalized(s));
		s="(S[+80.0]S)SSR";
		assertFalse(AmbiguousPeptideModSeq.isLocalized(s));
		s="(S[+80.0]S[+80.0])SSR";
		assertTrue(AmbiguousPeptideModSeq.isLocalized(s));
		s="(S[+80.0]S[+80.0]S)SR";
		assertFalse(AmbiguousPeptideModSeq.isLocalized(s));
		s="(S[+80.0])SS(S[+80.0])R";
		assertTrue(AmbiguousPeptideModSeq.isLocalized(s));
		s="(S[+80.0])S(SS[+80.0])R";
		assertFalse(AmbiguousPeptideModSeq.isLocalized(s));
	}
}
