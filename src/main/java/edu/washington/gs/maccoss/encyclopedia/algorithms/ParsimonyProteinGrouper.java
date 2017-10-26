package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroup;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ProteinGroupInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class ParsimonyProteinGrouper {
	public static Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>> groupProteins(ArrayList<PercolatorPeptide> passingPeptides, ArrayList<PercolatorPeptide> passingDecoys, float fdrThreshold) {
		ArrayList<PercolatorPeptide> allPeptides=new ArrayList<>();
		allPeptides.addAll(passingPeptides);
		allPeptides.addAll(passingDecoys);
		
		ArrayList<ProteinGroupInterface> allProteins=groupProteins(allPeptides);
		Collections.sort(allProteins);

		// calculate decoy/target FDR
		float[] fdrs=new float[allProteins.size()];
		int numberOfTargets=0;
		int numberOfDecoys=0;
		for (int i=allProteins.size()-1; i>=0; i--) {
			if (allProteins.get(i).isDecoy()) {
				numberOfDecoys++;
			} else {
				numberOfTargets++;
			}
			fdrs[i]=numberOfDecoys/(float)numberOfTargets;
		}
		
		// estimate q-values
		float minFDR=fdrs[0];
		for (int i=0; i<fdrs.length; i++) {
			if (fdrs[i]<minFDR) {
				minFDR=fdrs[i];
			} else {
				fdrs[i]=minFDR;
			}
		}
		
		ArrayList<PercolatorProteinGroup> keptTargets=new ArrayList<PercolatorProteinGroup>();
		ArrayList<PercolatorProteinGroup> keptDecoys=new ArrayList<PercolatorProteinGroup>();
		for (int i=allProteins.size()-1; i>=0; i--) {
			if (fdrs[i]<=fdrThreshold) {
				ProteinGroupInterface group=allProteins.get(i);
				if (group.isDecoy()) {
					keptDecoys.add(new PercolatorProteinGroup(group.getEquivalentAccessions(), group.getSequences(), fdrs[i], group.getPosteriorErrorProb()));
				} else {
					keptTargets.add(new PercolatorProteinGroup(group.getEquivalentAccessions(), group.getSequences(), fdrs[i], group.getPosteriorErrorProb()));
				}
			}
		}
		return new Pair<ArrayList<PercolatorProteinGroup>, ArrayList<PercolatorProteinGroup>>(keptTargets, keptDecoys);
	}
	
	public static ArrayList<ProteinGroupInterface> groupProteins(ArrayList<PercolatorPeptide> passingPeptides) {
		HashMap<String, Peptide> peptides=new HashMap<String, ParsimonyProteinGrouper.Peptide>();
		HashMap<String, Protein> proteins=new HashMap<String, ParsimonyProteinGrouper.Protein>();
		
		for (PercolatorPeptide percolatorPeptide : passingPeptides) {
			String sequence=PeptideUtils.getCorrectedMasses(PercolatorPeptide.getPeptideSequence(percolatorPeptide.getPsmID()));
			HashSet<String> accessions=PSMData.stringToAccessions(percolatorPeptide.getProteinIDs());
			
			Peptide peptide=peptides.get(sequence);
			if (peptide==null) {
				peptide=new Peptide(sequence, percolatorPeptide.getPosteriorErrorProb());
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

		ArrayList<Protein> sortedProteins=new ArrayList<ParsimonyProteinGrouper.Protein>(proteins.values());
		for (Protein protein : sortedProteins) {
			protein.recalculateNSP();
		}
		
		ArrayList<ProteinGroupInterface> keptProteins=new ArrayList<ProteinGroupInterface>();
		while (sortedProteins.size()>0) {
			Collections.sort(sortedProteins);
			Protein highestRankedProtein=sortedProteins.remove(sortedProteins.size()-1);
			if (highestRankedProtein.getNSP()==0.0f) {
				break;
			}
			ArrayList<Protein> equivalentProteins=highestRankedProtein.claimAllPeptides();
			
			HashSet<String> equivalentAccessions=new HashSet<String>();
			for (Protein protein : equivalentProteins) {
				equivalentAccessions.add(protein.accession);
				sortedProteins.remove(protein);
			}
			ArrayList<String> sequences=new ArrayList<String>();
			for (Peptide peptide : highestRankedProtein.peptides) {
				sequences.add(peptide.sequence);
			}
			keptProteins.add(new ProteinGroup(highestRankedProtein.getNSP(), highestRankedProtein.getMinPosterorErrorProbability(), new ArrayList<String>(equivalentAccessions), sequences));
		}
		return keptProteins;
	}

	static class Peptide {
		private final String sequence;
		private final float posteriorErrorProbability;
		private final ArrayList<Protein> proteins;

		public Peptide(String sequence, float posteriorErrorProbability) {
			this.sequence=sequence;
			this.posteriorErrorProbability=posteriorErrorProbability;
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
		private float minPosterorErrorProbability;

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
		
		public float getMinPosterorErrorProbability() {
			return minPosterorErrorProbability;
		}

		public void recalculateNSP() {
			nsp=0.0f;
			minPosterorErrorProbability=1.0f;
			for (Peptide peptide : peptides) {
				nsp+=(1.0f-peptide.posteriorErrorProbability);
				if (peptide.posteriorErrorProbability<minPosterorErrorProbability) {
					minPosterorErrorProbability=peptide.posteriorErrorProbability;
				}
			}
		}
		
		/**
		 * 
		 * @return returns identical proteins that contain the same peptides
		 */
		public ArrayList<Protein> claimAllPeptides() {
			HashSet<Protein> similarProteins=new HashSet<>();
			for (Peptide peptide : peptides) {
				similarProteins.addAll(peptide.proteins);
			}
			
			ArrayList<Protein> identicalProteins=new ArrayList<ParsimonyProteinGrouper.Protein>();
			for (Protein protein : similarProteins) {
				if (this.hasEquivalentPeptides(protein)) {
					identicalProteins.add(protein);
				}
			}
			
			for (Peptide peptide : peptides) {
				peptide.claimPeptide(this);
			}
			return identicalProteins;
		}
		
		public boolean hasEquivalentPeptides(Protein p) {
			if (peptides.size()!=p.peptides.size()) return false;
			
			HashSet<Peptide> sequences=new HashSet<>(peptides);
			for (Peptide peptide : p.peptides) {
				if (!sequences.contains(peptide)) return false;
			}
			return true;
		}

		@Override
		public int compareTo(Protein o) {
			if (o==null) return 1;
			int c=-Float.compare(minPosterorErrorProbability, o.minPosterorErrorProbability);
			if (c!=0) return c;
			c=Float.compare(nsp, o.nsp);
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
			return accession.equals(((Protein)obj).accession);
		}
	}
}
