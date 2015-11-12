package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.map.hash.TDoubleIntHashMap;

import java.util.ArrayList;
import java.util.Arrays;

public class BackgroundGenerator {
	/**
	 * TODO assumes no modifications!
	 * @param binBoundaries
	 * @param fasta
	 * @param params
	 */
	public static TDoubleIntHashMap[] generateBackground(double[] binBoundaries, ArrayList<FastaEntry> entries, SearchParameters params) {
		TDoubleIntHashMap[] binCounters=new TDoubleIntHashMap[binBoundaries.length-1];
		for (int i=0; i<binCounters.length; i++) {
			binCounters[i]=new TDoubleIntHashMap();
		}

		for (FastaEntry entry : entries) {
			ArrayList<String> peptides=params.getEnzyme().digestProtein(entry.getSequence(), params.getMinPeptideLength(), params.getMaxPeptideLength(), params.getMaxMissedCleavages());
			for (String sequence : peptides) {
				FragmentationModel model=new FragmentationModel(sequence);
				double[] ions=model.getPrimaryIons(params.getFragType());
				for (byte charge=params.getMinCharge(); charge<=params.getMaxCharge(); charge++) {
					double parentMZ=MassConstants.getChargedMass(sequence, charge);
					int index=Arrays.binarySearch(binBoundaries, parentMZ);
					if (index>=0) {
						// increment the lower index
						if (index>=binCounters.length) {
							continue;
						}
					} else {
						// insertion point
						index=-(index+1);
						if (index<=0||index>binCounters.length) {
							continue;
						}
						index--; // increment the lower index
					}
					
					for (double ion : ions) {
						binCounters[index].adjustOrPutValue(ion, 1, 1);
					}
				}
			}
		}
		
		return binCounters;
	}

}
