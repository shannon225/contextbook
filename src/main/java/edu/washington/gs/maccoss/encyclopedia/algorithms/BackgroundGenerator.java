package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptideDatabase;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.map.hash.TDoubleIntHashMap;

public class BackgroundGenerator {
	/**
	 * TODO assumes no modifications!
	 * 
	 * Note, does not use decoy proteins, which causes decoys to score somewhat
	 * higher due to having a slightly higher uniqueness factor. In a sense, the
	 * background just represents prior frequency assumptions based on nature,
	 * and it's actually better that decoys aren't considered.
	 * 
	 * @param binBoundaries
	 * @param fasta
	 * @param params
	 */

	public static Triplet<TDoubleIntHashMap[], ArrayList<String>[], HashSet<String>[]> generateBackground(double[] binBoundaries, boolean[] useBin, PeptideDatabase targets,
			HashSet<String> backgroundProteome, PecanSearchParameters parameters) {
		@SuppressWarnings("unchecked")
		HashSet<String>[] backgroundDecoys=new HashSet[binBoundaries.length];
		for (int i=0; i<backgroundDecoys.length; i++) {
			backgroundDecoys[i]=new HashSet<String>();
		}
		ArrayList<String> backgroundList=new ArrayList<String>(backgroundProteome);
		RandomGenerator.shuffle(backgroundList, 16807);
		for (String peptide : backgroundList) {
			for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
				double mz=parameters.getAAConstants().getChargedMass(peptide, charge);
				int index=Arrays.binarySearch(binBoundaries, mz);
				index=(-(index+1))-1;
				if (index>=0&&index<useBin.length) {
					if (useBin[index]) {
						if (backgroundDecoys[index].size()>2000) continue;

						String random=PeptideUtils.getDecoy(peptide, backgroundProteome, parameters);
						backgroundDecoys[index].add(random);
					}
				}
			}
		}

		if (parameters.isAddDecoysToBackgound()) {
			// add reverse targets to proteome
			for (FastaPeptideEntry entry : targets) {
				for (byte charge=parameters.getMinCharge(); charge<=parameters.getMaxCharge(); charge++) {
					String decoy=PeptideUtils.getSmartDecoy(entry.getSequence(), charge, backgroundProteome, parameters);
					backgroundProteome.add(decoy);
				}
			}
		}

		Pair<TDoubleIntHashMap[], ArrayList<String>[]> background=BackgroundGenerator.generateBackground(binBoundaries, useBin, backgroundProteome, backgroundDecoys, parameters);
		return new Triplet<TDoubleIntHashMap[], ArrayList<String>[], HashSet<String>[]>(background.x, background.y, backgroundDecoys);
	}

	public static Pair<TDoubleIntHashMap[], ArrayList<String>[]> generateBackground(double[] binBoundaries, boolean[] useBins, Collection<String> peptides, HashSet<String>[] backgroundDecoys,
			PecanSearchParameters params) {
		TDoubleIntHashMap[] binCounters=new TDoubleIntHashMap[binBoundaries.length-1];
		@SuppressWarnings("unchecked")
		HashSet<String>[] allPeptides=new HashSet[binBoundaries.length-1];
		for (int i=0; i<allPeptides.length; i++) {
			binCounters[i]=new TDoubleIntHashMap();
			allPeptides[i]=new HashSet<String>();
		}
		
		if (params.isAddDecoysToBackgound()) {
			// add bin-specific background
			for (int i=0; i<binCounters.length; i++) {
				allPeptides[i].addAll(backgroundDecoys[i]);
				for (String sequence : backgroundDecoys[i]) {
					for (byte charge=params.getMinCharge(); charge<=params.getMaxCharge(); charge++) {
						double parentMZ=params.getAAConstants().getChargedMass(sequence, charge);

						int index=Arrays.binarySearch(binBoundaries, parentMZ);
						index=(-(index+1))-1;
						if (index==i) {
							FragmentationModel model=new FragmentationModel(sequence, params.getAAConstants());
							double[] ions=model.getPrimaryIons(params.getFragType(), charge);
							for (double ion : ions) {
								binCounters[index].adjustOrPutValue(ion, 1, 1);
							}
						}
					}
				}
			}
		}
		
		// add actual proteome
		for (String sequence : peptides) {
			for (byte charge=params.getMinCharge(); charge<=params.getMaxCharge(); charge++) {
				double parentMZ=params.getAAConstants().getChargedMass(sequence, charge);

				int index=Arrays.binarySearch(binBoundaries, parentMZ);
				index=(-(index+1))-1;

				if (index<0||index>=binCounters.length||!useBins[index]) continue;

				FragmentationModel model=new FragmentationModel(sequence, params.getAAConstants());
				double[] ions=model.getPrimaryIons(params.getFragType(), charge);
				for (double ion : ions) {
					binCounters[index].adjustOrPutValue(ion, 1, 1);
				}
				allPeptides[index].add(sequence);
			}
		}

		@SuppressWarnings("unchecked")
		ArrayList<String>[] peptideArrays=new ArrayList[allPeptides.length];
		for (int i=0; i<peptideArrays.length; i++) {
			peptideArrays[i]=new ArrayList<String>(allPeptides[i]);
		}

		return new Pair<TDoubleIntHashMap[], ArrayList<String>[]>(binCounters, peptideArrays);
	}

	/**
	 * TODO assumes no modifications!
	 * 
	 * Note, does not use decoy proteins, which causes decoys to score somewhat
	 * higher due to having a slightly higher uniqueness factor. In a sense, the
	 * background just represents prior frequency assumptions based on nature,
	 * and it's actually better that decoys aren't considered.
	 * 
	 * @param binBoundaries
	 * @param fasta
	 * @param params
	 */
	public static Pair<TDoubleIntHashMap[], ArrayList<String>[]> generateBackground(double[] binBoundaries, Collection<FastaEntryInterface> entries, boolean digest, PecanSearchParameters params) {
		@SuppressWarnings("unchecked")
		HashSet<String>[] backgroundDecoys=new HashSet[binBoundaries.length];
		boolean[] useBins=new boolean[binBoundaries.length];
		for (int i=0; i<backgroundDecoys.length; i++) {
			backgroundDecoys[i]=new HashSet<String>();
			useBins[i]=true;
		}

		HashSet<String> sequences=new HashSet<String>();
		for (FastaEntryInterface entry : entries) {
			ArrayList<String> peptides;
			if (digest) {
				peptides=params.getEnzyme().digestProtein(entry.getSequence(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxMissedCleavages(), params.getAAConstants());
			} else {
				peptides=new ArrayList<String>();
				peptides.add(entry.getSequence());
			}
			for (String sequence : peptides) {
				sequences.add(sequence);
			}
		}
		return generateBackground(binBoundaries, useBins, sequences, backgroundDecoys, params);
	}

}
