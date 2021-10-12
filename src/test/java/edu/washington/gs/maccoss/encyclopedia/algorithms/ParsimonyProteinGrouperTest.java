package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PercolatorReader;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.TestCase;

public class ParsimonyProteinGrouperTest extends TestCase {
	public static void main(String[] args) {
		File outputFile=new File("/Volumes/BriansSSD/hela_serum_timecourse/hela_serum_timecourse_wide_window_combined_concatenated_results.txt");
		File decoyFile=new File("/Volumes/BriansSSD/hela_serum_timecourse/hela_serum_timecourse_wide_window_combined_concatenated_decoy.txt");
		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());
		ArrayList<PercolatorPeptide> targets=PercolatorReader.getPassingPeptidesFromTSV(outputFile, 0.01f, aaConstants, true).x;
		ArrayList<PercolatorPeptide> decoys=PercolatorReader.getPassingPeptidesFromTSV(decoyFile, 0.01f, aaConstants, true).x;
		
		System.out.println("NONE:    "+ParsimonyProteinGrouper.groupPercolatorProteins(targets, aaConstants).size());
		System.out.println("10%FDR:  "+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.1f, aaConstants).x.size());
		System.out.println(" 5%FDR:  "+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.05f, aaConstants).x.size());
		System.out.println(" 1%FDR:  "+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.01f, aaConstants).x.size());
		System.out.println("0.1%FDR: "+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.001f, aaConstants).x.size());
		System.out.println("0.01%FDR:"+ParsimonyProteinGrouper.groupProteins(targets, decoys, 0.0001f, aaConstants).x.size());
	}

	public void testProteinGrouper() {
		ArrayList<PercolatorPeptide> peptides=new ArrayList<PercolatorPeptide>();
		final AminoAcidConstants aaConstants = new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap());

		int randomPeptides=10;
		for (int i=0; i<randomPeptides; i++) {
			String psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, Optional.empty(), false, RandomGenerator.randomSequence(i), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+RandomGenerator.randomSequence(i), 0.0f, 0.0f, aaConstants));
		}
		ArrayList<ProteinGroupInterface> proteins=ParsimonyProteinGrouper.groupPercolatorProteins(peptides, aaConstants);
		assertEquals(randomPeptides, proteins.size());
		
		peptides.clear();
		for (int i=0; i<randomPeptides; i++) {
			String accession=RandomGenerator.randomSequence(i);
			String psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, Optional.empty(), false, RandomGenerator.randomSequence(i), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+accession, 0.0f, 0.0f, aaConstants));
			psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, Optional.empty(), false, RandomGenerator.randomSequence(i+10), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+accession, 0.0f, 0.0f, aaConstants));
		}
		proteins=ParsimonyProteinGrouper.groupPercolatorProteins(peptides, aaConstants);
		assertEquals(randomPeptides, proteins.size());
		
		for (ProteinGroupInterface proteinGroup : proteins) {
			assertEquals(1, proteinGroup.getEquivalentAccessions().size());
		}
		
		peptides.clear();
		for (int i=0; i<randomPeptides; i++) {
			String accession=RandomGenerator.randomSequence(i);
			String altAccession=RandomGenerator.randomSequence(i)+"-alt";
			String psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, Optional.empty(), false, RandomGenerator.randomSequence(i), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+accession, 0.0f, 0.0f, aaConstants));
			peptides.add(new PercolatorPeptide(psmID, ">"+altAccession, 0.0f, 0.0f, aaConstants));
			psmID=PercolatorPeptide.getPSMID("FILE", 0.0f, Optional.empty(), false, RandomGenerator.randomSequence(i+10), (byte)2);
			peptides.add(new PercolatorPeptide(psmID, ">"+accession, 0.0f, 0.0f, aaConstants));
			peptides.add(new PercolatorPeptide(psmID, ">"+altAccession, 0.0f, 0.0f, aaConstants));
		}
		proteins=ParsimonyProteinGrouper.groupPercolatorProteins(peptides, aaConstants);
		assertEquals(randomPeptides, proteins.size());
		
		for (ProteinGroupInterface proteinGroup : proteins) {
			assertEquals(2, proteinGroup.getEquivalentAccessions().size());
		}
	}
}
