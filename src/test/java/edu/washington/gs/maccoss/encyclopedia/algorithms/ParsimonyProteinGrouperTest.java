package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import junit.framework.TestCase;

public class ParsimonyProteinGrouperTest extends TestCase {
	public static void main(String[] args) {
		File outputFile=new File("/Volumes/BriansSSD/pecan/group_concatenated_results.txt");
		File decoyFile=new File("/Volumes/BriansSSD/pecan/group_concatenated_decoy.txt");
		ArrayList<PercolatorPeptide> targets=PercolatorReader.getPassingPeptidesFromTSV(outputFile, 0.01f, true).x;
		ArrayList<PercolatorPeptide> decoys=PercolatorReader.getPassingPeptidesFromTSV(decoyFile, 0.01f, true).x;
		
		System.out.println("NONE:    "+ParsimonyProteinGrouper.groupProteins(targets).size());
		System.out.println("10%FDR:  "+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.1f).x.size());
		System.out.println(" 5%FDR:  "+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.05f).x.size());
		System.out.println(" 1%FDR:  "+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.01f).x.size());
		System.out.println("0.1%FDR: "+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.001f).x.size());
		System.out.println("0.01%FDR:"+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.0001f).x.size());
	}
	
	public void testProteinGrouper() {
		ArrayList<PercolatorPeptide> peptides=new ArrayList<PercolatorPeptide>();
		
		int randomPeptides=10;
		for (int i=0; i<randomPeptides; i++) {
			String psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, false, RandomGenerator.randomSequence(i), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+RandomGenerator.randomSequence(i), 0.0f, 0.0f));
		}
		ArrayList<ProteinGroupInterface> proteins=ParsimonyProteinGrouper.groupProteins(peptides);
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
		
		for (ProteinGroupInterface proteinGroup : proteins) {
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
		
		for (ProteinGroupInterface proteinGroup : proteins) {
			assertEquals(2, proteinGroup.getEquivalentAccessions().size());
		}
	}
}
