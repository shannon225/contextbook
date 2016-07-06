package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;

public class ParsimonyProteinGrouper {
	public static void main(String[] args) throws Exception {
		/*File f=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.mzml.encyclopedia.txt");
		f=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/121115_BCS_HeLa_24mz_400_1000.dia.encyclopedia.txt");
		f=new File("/Users/searleb/Documents/projects/encyclopedia/mzml/deep_hela/concatenated_results.txt");
		
		Logger.logLine("Starting reading Percolator result...");
		ArrayList<PercolatorPeptide> passingPeptidesFromTSV=PercolatorReader.getPassingPeptidesFromTSV(f, 0.01f);
		Logger.logLine("Starting grouping proteins...");
		groupProtein(passingPeptidesFromTSV);
		Logger.logLine("Finished!");*/
		
		File libraryFile=new File("/Users/searleb/Documents/projects/encyclopedia/HumanTotalProteome/HeLa.elib");
		libraryFile=new File("/Users/searleb/Documents/school/projects/may_asms/yeast/YeastProteome.elib");
		libraryFile=new File("/Users/searleb/Documents/projects/phosphopedia/VillenJ_Exactive_HumanPhosphoproteome.elib");
		LibraryFile file=new LibraryFile();
		file.openFile(libraryFile);
		
		ArrayList<PercolatorPeptide> peptides=new ArrayList<PercolatorPeptide>();
		for (LibraryEntry entry : file.getAllEntries(false)) {
			peptides.add(entry.getPSMData());
		}
		Logger.logLine("Starting grouping proteins...");
		groupProtein(peptides);
		Logger.logLine("Finished!");
	}
	
	public static ArrayList<ScoredObject<String>> groupProtein(ArrayList<PercolatorPeptide> passingPeptides) {
		HashMap<String, Peptide> peptides=new HashMap<String, ParsimonyProteinGrouper.Peptide>();
		HashMap<String, Protein> proteins=new HashMap<String, ParsimonyProteinGrouper.Protein>();
		
		for (PercolatorPeptide percolatorPeptide : passingPeptides) {
			String sequence=PercolatorPeptide.getPeptideSequence(percolatorPeptide.getPsmID());
			HashSet<String> accessions=PSMData.stringToAccessions(percolatorPeptide.getProteinIDs());
			
			Peptide peptide=peptides.get(sequence);
			if (peptide==null) {
				peptide=new Peptide(sequence, 1.0f-percolatorPeptide.getPosteriorErrorProb());
				peptides.put(sequence, peptide);
			}
			
			for (String accession : accessions) {
				Protein protein=proteins.get(accession);
				if (protein==null) {
					protein=new Protein(accession);
					proteins.put(accession, protein);
				}
				peptide.addProtein(protein);
				protein.addPeptide(peptide);
			}
		}

		Logger.logLine(proteins.size()+" total accessions from "+peptides.size()+" peptides...");

		ArrayList<Protein> sortedProteins=new ArrayList<ParsimonyProteinGrouper.Protein>(proteins.values());
		for (Protein protein : sortedProteins) {
			protein.recalcualteNSP();
		}
		ArrayList<ScoredObject<String>> keptProteins=new ArrayList<ScoredObject<String>>();
		while (sortedProteins.size()>0) {
			Collections.sort(sortedProteins);
			Protein highestRankedProtein=sortedProteins.remove(sortedProteins.size()-1);
			if (highestRankedProtein.getNSP()==0.0f) {
				break;
			}
			keptProteins.add(new ScoredObject<String>(highestRankedProtein.getNSP(), highestRankedProtein.accession));
			highestRankedProtein.claimAllPeptides();
		}
		
		Logger.logLine(keptProteins.size()+" parsimonious proteins from "+peptides.size()+" peptides");
		return keptProteins;
	}

	static class Peptide {
		private final String sequence;
		private final float probability;
		private final ArrayList<Protein> proteins;

		public Peptide(String sequence, float probability) {
			this.sequence=sequence;
			this.probability=probability;
			this.proteins=new ArrayList<ParsimonyProteinGrouper.Protein>();
		}

		public void addProtein(Protein protein) {
			proteins.add(protein);
		}

		public void claimPeptide(Protein claimer) {
			for (Protein protein : proteins) {
				if (protein!=claimer) {
					protein.removePeptide(this);
					protein.recalcualteNSP();
				}
			}
		}

		@Override
		public int hashCode() {
			return sequence.hashCode();
		}

		@Override
		public boolean equals(Object obj) {
			if (obj==null) return false;
			return sequence.equals(((Peptide)obj).sequence);
		}
	}

	static class Protein implements Comparable<Protein> {
		private final String accession;
		private final ArrayList<Peptide> peptides;
		private float nsp;

		public Protein(String accession) {
			this.accession=accession;
			this.peptides=new ArrayList<ParsimonyProteinGrouper.Peptide>();
		}

		public void addPeptide(Peptide peptide) {
			peptides.add(peptide);
		}

		public void removePeptide(Peptide peptide) {
			peptides.remove(peptide);
		}

		public float getNSP() {
			return nsp;
		}

		public void recalcualteNSP() {
			nsp=0.0f;
			for (Peptide peptide : peptides) {
				nsp+=peptide.probability;
			}
		}
		
		public void claimAllPeptides() {
			for (Peptide peptide : peptides) {
				peptide.claimPeptide(this);
			}
		}

		@Override
		public int compareTo(Protein o) {
			if (o==null) return 1;
			int c=Float.compare(nsp, o.nsp);
			if (c!=0) return c;
			return accession.compareTo(o.accession);
		}

		@Override
		public int hashCode() {
			return accession.hashCode();
		}

		@Override
		public boolean equals(Object obj) {
			if (obj==null) return false;
			return accession.equals(((Peptide)obj).sequence);
		}
	}
}
