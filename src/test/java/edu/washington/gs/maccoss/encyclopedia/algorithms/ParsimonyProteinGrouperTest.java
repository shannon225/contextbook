package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import junit.framework.TestCase;

public class ParsimonyProteinGrouperTest extends TestCase {
	public void testProteinGrouper() {
		ArrayList<PercolatorPeptide> peptides=new ArrayList<PercolatorPeptide>();
		
		int randomPeptides=10;
		for (int i=0; i<randomPeptides; i++) {
			String psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, false, RandomGenerator.randomSequence(i), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+RandomGenerator.randomSequence(i), 0.0f, 0.0f));
		}
		ArrayList<ProteinGroup> proteins=ParsimonyProteinGrouper.groupProteins(peptides);
		assertEquals(randomPeptides, proteins.size());
		
		peptides.clear();
		for (int i=0; i<randomPeptides; i++) {
			String accession=RandomGenerator.randomSequence(i);
			String psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, false, RandomGenerator.randomSequence(i), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+accession, 0.0f, 0.0f));
			psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, false, RandomGenerator.randomSequence(i+10), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+accession, 0.0f, 0.0f));
		}
		proteins=ParsimonyProteinGrouper.groupProteins(peptides);
		assertEquals(randomPeptides, proteins.size());
		
		for (ProteinGroup proteinGroup : proteins) {
			assertEquals(1, proteinGroup.getEquivalentAccessions().size());
		}
		
		peptides.clear();
		for (int i=0; i<randomPeptides; i++) {
			String accession=RandomGenerator.randomSequence(i);
			String altAccession=RandomGenerator.randomSequence(i)+"-alt";
			String psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, false, RandomGenerator.randomSequence(i), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+accession, 0.0f, 0.0f));
			peptides.add(new PercolatorPeptide(psmID, ">"+altAccession, 0.0f, 0.0f));
			psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, false, RandomGenerator.randomSequence(i+10), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+accession, 0.0f, 0.0f));
			peptides.add(new PercolatorPeptide(psmID, ">"+altAccession, 0.0f, 0.0f));
		}
		proteins=ParsimonyProteinGrouper.groupProteins(peptides);
		assertEquals(randomPeptides, proteins.size());
		
		for (ProteinGroup proteinGroup : proteins) {
			assertEquals(2, proteinGroup.getEquivalentAccessions().size());
		}
	}
}
