package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

public class AmbiguousPeptideModSeqTest extends TestCase {
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"));
	public static void main(String[] args) {
		System.out.println(AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.96633]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.96633]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getLeftAmbiguity("SS[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getLeftAmbiguity("SSS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getRightAmbiguity("SSS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		
	}
	
	public void testAmbiguity() {;
		assertEquals("<S[+79.96633])SSR", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.96633]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]SS>R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.96633]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("<SS[+79.96633])SR", AmbiguousPeptideModSeq.getLeftAmbiguity("SS[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.96633]S>R", AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("<SSS[+79.96633])R", AmbiguousPeptideModSeq.getLeftAmbiguity("SSS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("SS(S[+79.96633]>R", AmbiguousPeptideModSeq.getRightAmbiguity("SSS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		

		assertEquals("<S[+79.96633])<S[+79.96633])SR", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.96633]S[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]>(S[+79.96633]S>R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.96633]S[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("<SS[+79.96633])<S[+79.96633])R", AmbiguousPeptideModSeq.getLeftAmbiguity("SS[+79.96633]S[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.96633]>(S[+79.96633]>R", AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.96633]S[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("<S[+79.96633])<SS[+79.96633])R", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.96633]SS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]S>(S[+79.96633]>R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.96633]SS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());

		assertEquals("(S[+79.96633])SSR", AmbiguousPeptideModSeq.getUnambigous("S[+79.96633]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.96633])SR", AmbiguousPeptideModSeq.getUnambigous("SS[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("SS(S[+79.96633])R", AmbiguousPeptideModSeq.getUnambigous("SSS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633])(S[+79.96633])SR", AmbiguousPeptideModSeq.getUnambigous("S[+79.96633]S[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.96633])(S[+79.96633])R", AmbiguousPeptideModSeq.getUnambigous("SS[+79.96633]S[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633])S(S[+79.96633])R", AmbiguousPeptideModSeq.getUnambigous("S[+79.96633]SS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());

		assertEquals("(S[+79.96633]SS)R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.96633]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SS[+79.96633]S)R", AmbiguousPeptideModSeq.getFullyAmbiguous("SS[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SSS[+79.96633])R", AmbiguousPeptideModSeq.getFullyAmbiguous("SSS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]S[+79.96633]S)R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.96633]S[+79.96633]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SS[+79.96633]S[+79.96633])R", AmbiguousPeptideModSeq.getFullyAmbiguous("SS[+79.96633]S[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.96633]SS[+79.96633])R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.96633]SS[+79.96633]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
	}
	
	public void testAmbiguousSorting() {
		String[] peptides=new String[] {
				"SS[+79.96633]SR",
				"S[+79.96633]S[+79.96633]SR",
				"S[+79.96633]S[+79.96633]S[+79.96633]R",
				"S[+79.96633]QWEITS[+79.96633]GLKDSS[+79.96633]R",
				"SQWEIT[+79.96633]SGLKDS[+79.96633]S[+79.96633]R"
		};
		int[] diffs=new int[] {
				0, -1, 0, -8, 9
		};
		for (int i=0; i<peptides.length; i++) {
			String targetPeptide=peptides[i];
			AmbiguousPeptideModSeq left=AmbiguousPeptideModSeq.getLeftAmbiguity(targetPeptide, PeptideModification.phosphorylation,PARAMETERS.getAAConstants());
			AmbiguousPeptideModSeq right=AmbiguousPeptideModSeq.getRightAmbiguity(targetPeptide, PeptideModification.phosphorylation,PARAMETERS.getAAConstants());
			assertEquals(diffs[i], left.numAmbigousResidues()-right.numAmbigousResidues());
		}
	}
	
	public void testRemove() {
		String targetPeptide="S[+79.96633]SSR";

		String[] unambiguousPeptides=new String[] {
				"S[+79.96633]SSR",
				"SS[+79.96633]SR",
				"SSS[+79.96633]R"
		};
		String[] expectedResult=new String[] {
				"S(S[+79.96633]S)R", // moved!
				"(S[+79.96633]SS)R",
				"(S[+79.96633]S)SR"
		};
		for (int i=0; i<unambiguousPeptides.length; i++) {
			AmbiguousPeptideModSeq s=AmbiguousPeptideModSeq.getFullyAmbiguous(targetPeptide, PeptideModification.phosphorylation, PARAMETERS.getAAConstants());
			AmbiguousPeptideModSeq unambiguous=AmbiguousPeptideModSeq.getUnambigous(unambiguousPeptides[i], PeptideModification.phosphorylation, PARAMETERS.getAAConstants());
			
			assertEquals(expectedResult[i], s.removeAmbiguity(unambiguous).get().getPeptideAnnotation());
		}
	}
	
	public void testSets() {
		String[] peptides=new String[] {
				"SS[+79.96633]SR",
				"S[+79.96633]S[+79.96633]SR",
				"S[+79.96633]S[+79.96633]S[+79.96633]R",
				"S[+79.96633]QWEITS[+79.96633]GLKDSS[+79.96633]R"
		};
		for (String targetPeptide : peptides) {
			AmbiguousPeptideModSeq s=AmbiguousPeptideModSeq.getRightAmbiguity(targetPeptide, PeptideModification.phosphorylation,PARAMETERS.getAAConstants());
			TIntHashSet[] sets=s.getAmbiguityGroups();
			assertEquals(PeptideUtils.getNumberOfMods(targetPeptide, 80), sets.length);
			int[] setData=AmbiguousPeptideModSeq.getModificationGroupsFromSets(sets, s.length());
			int[] groups=s.getModificationGroup();
			assertEquals(groups.length, setData.length);
			for (int i=0; i<groups.length; i++) {
				assertEquals(groups[i], setData[i]);
			}
		}
	}
	
	public void testIsLocalized() {
		String s="(S[+80.0])SSSR";
		assertTrue(AmbiguousPeptideModSeq.isLocalized(s, PeptideModification.phosphorylation));
		s="(S[+80.0]S)SSR";
		assertFalse(AmbiguousPeptideModSeq.isLocalized(s, PeptideModification.phosphorylation));
		s="(S[+80.0]S[+80.0])SSR";
		assertTrue(AmbiguousPeptideModSeq.isLocalized(s, PeptideModification.phosphorylation));
		s="(S[+80.0]S[+80.0]S)SR";
		assertFalse(AmbiguousPeptideModSeq.isLocalized(s, PeptideModification.phosphorylation));
		s="(S[+80.0])SS(S[+80.0])R";
		assertTrue(AmbiguousPeptideModSeq.isLocalized(s, PeptideModification.phosphorylation));
		s="(S[+80.0])S(SS[+80.0])R";
		assertFalse(AmbiguousPeptideModSeq.isLocalized(s, PeptideModification.phosphorylation));
	}
	
	public void testIsLocalizedAtEnd() {
		assertTrue(AmbiguousPeptideModSeq.isLocalizedAtEnd("<S[+80.0])SSSR", PeptideModification.phosphorylation));
		assertTrue(AmbiguousPeptideModSeq.isLocalizedAtEnd("SSS(S[+80.0]>R", PeptideModification.phosphorylation));
		assertFalse(AmbiguousPeptideModSeq.isLocalizedAtEnd("(S[+80.0]>SSSR", PeptideModification.phosphorylation));
		assertFalse(AmbiguousPeptideModSeq.isLocalizedAtEnd("SSS<S[+80.0])R", PeptideModification.phosphorylation));
	}
}
