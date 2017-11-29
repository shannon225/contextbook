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
	private static final SearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"), false);
	public static void main(String[] args) {
		System.out.println(AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.966331]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.966331]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getLeftAmbiguity("SS[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getLeftAmbiguity("SSS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		System.out.println(AmbiguousPeptideModSeq.getRightAmbiguity("SSS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()));
		
	}
	
	public void testRetrieval() {
		String annotationString="KG<S[+79.966331])GDYMPMSPK";
		AmbiguousPeptideModSeq target=AmbiguousPeptideModSeq.getLeftAmbiguity("KGS[+79.966331]GDYMPMSPK", PeptideModification.phosphorylation, PARAMETERS.getAAConstants());
		assertEquals(annotationString, target.getPeptideAnnotation());
		
		assertRetrieval(annotationString);
		assertRetrieval("KG<S[+79.966331]GDY)MPMSPK");
		assertRetrieval("KG(S[+79.966331]GDY>MPMSPK");
		assertRetrieval("KG(SGDY[+79.966331]>MPMSPK");
		assertRetrieval("KGSGDYMPM(S[+79.966331]>PK");
	}

	private void assertRetrieval(String annotationString) {
		AminoAcidConstants aminoAcidConstants = new AminoAcidConstants();
		AmbiguousPeptideModSeq seq=AmbiguousPeptideModSeq.getAmbiguousPeptideModSeq(annotationString, PeptideModification.phosphorylation, aminoAcidConstants);
		assertEquals(annotationString, seq.getPeptideAnnotation());
	}

	
	public void testRemoveMore() {
		String peptideModSeq="KGS[+79.966331]GDYMPMSPK";
		AmbiguousPeptideModSeq right=AmbiguousPeptideModSeq.getRightAmbiguity(peptideModSeq, PeptideModification.phosphorylation,PARAMETERS.getAAConstants());
		assertEquals("KG(S[+79.966331]GDYMPMS>PK", right.getPeptideAnnotation());
		
		AmbiguousPeptideModSeq unambig=AmbiguousPeptideModSeq.getUnambigous("KGSGDYMPMS[+79.966331]PK", PeptideModification.phosphorylation,PARAMETERS.getAAConstants());
		AmbiguousPeptideModSeq unambig2=AmbiguousPeptideModSeq.getUnambigous("KGSGDY[+79.966331]MPMSPK", PeptideModification.phosphorylation,PARAMETERS.getAAConstants());
		
		assertEquals("KG(S[+79.966331]GDY>MPMSPK", right.removeAmbiguity(unambig).get().getPeptideAnnotation());
		assertEquals("KG(S[+79.966331]>GDYMPMSPK", right.removeAmbiguity(unambig, unambig2).get().getPeptideAnnotation());
	}
	
	public void testAmbiguity() {;
		assertEquals("<S[+79.966331])SSR", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.966331]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.966331]SS>R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.966331]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("<SS[+79.966331])SR", AmbiguousPeptideModSeq.getLeftAmbiguity("SS[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.966331]S>R", AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("<SSS[+79.966331])R", AmbiguousPeptideModSeq.getLeftAmbiguity("SSS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("SS(S[+79.966331]>R", AmbiguousPeptideModSeq.getRightAmbiguity("SSS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		

		assertEquals("<S[+79.966331])<S[+79.966331])SR", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.966331]S[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.966331]>(S[+79.966331]S>R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.966331]S[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("<SS[+79.966331])<S[+79.966331])R", AmbiguousPeptideModSeq.getLeftAmbiguity("SS[+79.966331]S[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.966331]>(S[+79.966331]>R", AmbiguousPeptideModSeq.getRightAmbiguity("SS[+79.966331]S[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("<S[+79.966331])<SS[+79.966331])R", AmbiguousPeptideModSeq.getLeftAmbiguity("S[+79.966331]SS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.966331]S>(S[+79.966331]>R", AmbiguousPeptideModSeq.getRightAmbiguity("S[+79.966331]SS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());

		assertEquals("(S[+79.966331])SSR", AmbiguousPeptideModSeq.getUnambigous("S[+79.966331]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.966331])SR", AmbiguousPeptideModSeq.getUnambigous("SS[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("SS(S[+79.966331])R", AmbiguousPeptideModSeq.getUnambigous("SSS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.966331])(S[+79.966331])SR", AmbiguousPeptideModSeq.getUnambigous("S[+79.966331]S[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("S(S[+79.966331])(S[+79.966331])R", AmbiguousPeptideModSeq.getUnambigous("SS[+79.966331]S[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.966331])S(S[+79.966331])R", AmbiguousPeptideModSeq.getUnambigous("S[+79.966331]SS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());

		assertEquals("(S[+79.966331]SS)R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.966331]SSR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SS[+79.966331]S)R", AmbiguousPeptideModSeq.getFullyAmbiguous("SS[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SSS[+79.966331])R", AmbiguousPeptideModSeq.getFullyAmbiguous("SSS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.966331]S[+79.966331]S)R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.966331]S[+79.966331]SR", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(SS[+79.966331]S[+79.966331])R", AmbiguousPeptideModSeq.getFullyAmbiguous("SS[+79.966331]S[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
		assertEquals("(S[+79.966331]SS[+79.966331])R", AmbiguousPeptideModSeq.getFullyAmbiguous("S[+79.966331]SS[+79.966331]R", PeptideModification.phosphorylation,PARAMETERS.getAAConstants()).getPeptideAnnotation());
	}
	
	public void testAmbiguousSorting() {
		String[] peptides=new String[] {
				"SS[+79.966331]SR",
				"S[+79.966331]S[+79.966331]SR",
				"S[+79.966331]S[+79.966331]S[+79.966331]R",
				"S[+79.966331]QWEITS[+79.966331]GLKDSS[+79.966331]R",
				"SQWEIT[+79.966331]SGLKDS[+79.966331]S[+79.966331]R"
		};
		int[] diffs=new int[] {
				0, -1, 0, -8, 9
		};
		for (int i=0; i<peptides.length; i++) {
			String targetPeptide=peptides[i];
			AmbiguousPeptideModSeq left=AmbiguousPeptideModSeq.getLeftAmbiguity(targetPeptide, PeptideModification.phosphorylation,PARAMETERS.getAAConstants());
			AmbiguousPeptideModSeq right=AmbiguousPeptideModSeq.getRightAmbiguity(targetPeptide, PeptideModification.phosphorylation,PARAMETERS.getAAConstants());
			assertEquals(diffs[i], left.getNumAmbigousResidues()-right.getNumAmbigousResidues());
		}
	}
	
	public void testRemove() {
		String targetPeptide="S[+79.966331]SSR";

		String[] unambiguousPeptides=new String[] {
				"S[+79.966331]SSR",
				"SS[+79.966331]SR",
				"SSS[+79.966331]R"
		};
		String[] expectedResult=new String[] {
				"S(S[+79.966331]S)R", // moved!
				"(S[+79.966331]SS)R",
				"(S[+79.966331]S)SR"
		};
		for (int i=0; i<unambiguousPeptides.length; i++) {
			AmbiguousPeptideModSeq s=AmbiguousPeptideModSeq.getFullyAmbiguous(targetPeptide, PeptideModification.phosphorylation, PARAMETERS.getAAConstants());
			AmbiguousPeptideModSeq unambiguous=AmbiguousPeptideModSeq.getUnambigous(unambiguousPeptides[i], PeptideModification.phosphorylation, PARAMETERS.getAAConstants());
			
			assertEquals(expectedResult[i], s.removeAmbiguity(unambiguous).get().getPeptideAnnotation());
		}
	}
	
	public void testSets() {
		String[] peptides=new String[] {
				"SS[+79.966331]SR",
				"S[+79.966331]S[+79.966331]SR",
				"S[+79.966331]S[+79.966331]S[+79.966331]R",
				"S[+79.966331]QWEITS[+79.966331]GLKDSS[+79.966331]R"
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
		assertTrue(AmbiguousPeptideModSeq.isSiteSpecific(s, PeptideModification.phosphorylation));
		s="(S[+80.0]S)SSR";
		assertFalse(AmbiguousPeptideModSeq.isSiteSpecific(s, PeptideModification.phosphorylation));
		s="(S[+80.0]S[+80.0])SSR";
		assertTrue(AmbiguousPeptideModSeq.isSiteSpecific(s, PeptideModification.phosphorylation));
		s="(S[+80.0]S[+80.0]S)SR";
		assertFalse(AmbiguousPeptideModSeq.isSiteSpecific(s, PeptideModification.phosphorylation));
		s="(S[+80.0])SS(S[+80.0])R";
		assertTrue(AmbiguousPeptideModSeq.isSiteSpecific(s, PeptideModification.phosphorylation));
		s="(S[+80.0])S(SS[+80.0])R";
		assertFalse(AmbiguousPeptideModSeq.isSiteSpecific(s, PeptideModification.phosphorylation));
	}
	
	public void testIsCompletelyAmbiguous() {
		String s="(S[+80.0])SSSR";
		assertFalse(AmbiguousPeptideModSeq.isCompletelyAmbiguous(s, PeptideModification.phosphorylation));
		s="(S[+80.0]SSS)R";
		assertTrue(AmbiguousPeptideModSeq.isCompletelyAmbiguous(s, PeptideModification.phosphorylation));
		s="(S[+80.0]S[+80.0])SSR";
		assertFalse(AmbiguousPeptideModSeq.isCompletelyAmbiguous(s, PeptideModification.phosphorylation));
		s="(S[+80.0]S[+80.0]S)SR";
		assertFalse(AmbiguousPeptideModSeq.isCompletelyAmbiguous(s, PeptideModification.phosphorylation));
		s="(S[+80.0]SSS[+80.0])R";
		assertTrue(AmbiguousPeptideModSeq.isCompletelyAmbiguous(s, PeptideModification.phosphorylation));
		s="(S[+80.0])S(SS[+80.0])R";
		assertFalse(AmbiguousPeptideModSeq.isCompletelyAmbiguous(s, PeptideModification.phosphorylation));
	}
	
	public void testIsLocalizedAtEnd() {
		assertTrue(AmbiguousPeptideModSeq.isSiteSpecificAtEnd("<S[+80.0])SSSR", PeptideModification.phosphorylation));
		assertTrue(AmbiguousPeptideModSeq.isSiteSpecificAtEnd("SSS(S[+80.0]>R", PeptideModification.phosphorylation));
		assertFalse(AmbiguousPeptideModSeq.isSiteSpecificAtEnd("(S[+80.0]>SSSR", PeptideModification.phosphorylation));
		assertFalse(AmbiguousPeptideModSeq.isSiteSpecificAtEnd("SSS<S[+80.0])R", PeptideModification.phosphorylation));
	}
}
