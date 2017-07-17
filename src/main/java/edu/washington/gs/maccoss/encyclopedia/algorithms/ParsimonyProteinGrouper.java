package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class ParsimonyProteinGrouper {
	public static ArrayList<ProteinGroup> groupProteins(ArrayList<PercolatorPeptide> passingPeptides) {
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
			protein.recalculateNSP();
		}
		ArrayList<ProteinGroup> keptProteins=new ArrayList<ProteinGroup>();
		while (sortedProteins.size()>0) {
			Collections.sort(sortedProteins);
			Protein highestRankedProtein=sortedProteins.remove(sortedProteins.size()-1);
			if (highestRankedProtein.getNSP()==0.0f) {
				break;
			}
			float nspScore=highestRankedProtein.getNSP();
			ArrayList<Protein> equivalentProteins=highestRankedProtein.claimAllPeptides();
			
			HashSet<String> equivalentAccessions=new HashSet<String>();
			for (Protein protein : equivalentProteins) {
				equivalentAccessions.add(protein.accession);
			}
			ArrayList<String> sequences=new ArrayList<String>();
			for (Peptide peptide : highestRankedProtein.peptides) {
				sequences.add(peptide.sequence);
			}
			keptProteins.add(new ProteinGroup(nspScore, new ArrayList<String>(equivalentAccessions), sequences));
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
					protein.recalculateNSP();
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

		public void recalculateNSP() {
			nsp=0.0f;
			for (Peptide peptide : peptides) {
				nsp+=peptide.probability;
			}
		}
		
		/**
		 * 
		 * @return returns identical proteins that contain the same peptides
		 */
		public ArrayList<Protein> claimAllPeptides() {
			ArrayList<Protein> identicalProteins=new ArrayList<ParsimonyProteinGrouper.Protein>();
			boolean first=true;
			for (Peptide peptide : peptides) {
				if (first) {
					identicalProteins.addAll(peptide.proteins);
				} else {
					ArrayList<Protein> toBeRemoved=new ArrayList<ParsimonyProteinGrouper.Protein>();
					for (Protein protein : identicalProteins) {
						if (!peptide.proteins.contains(protein)) {
							toBeRemoved.add(protein);
						}
					}
					identicalProteins.removeAll(toBeRemoved);
				}
				peptide.claimPeptide(this);
			}
			return identicalProteins;
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
