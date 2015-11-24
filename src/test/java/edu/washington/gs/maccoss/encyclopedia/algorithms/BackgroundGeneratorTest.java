package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import junit.framework.TestCase;

public class BackgroundGeneratorTest extends TestCase {
	private static final SearchParameters PARAMETERS=new SearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"));

	public void testGenerateBackground() {
		InputStream is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "contam_correctNL.fasta");
		TDoubleArrayList bins=new TDoubleArrayList();
		for (double i=400.0; i<=900.0; i+=5.0) {
			bins.add(i);
		}
		double[] binArray=bins.toArray();
		TDoubleIntHashMap[] binCounters=BackgroundGenerator.generateBackground(binArray, entries, PARAMETERS).x;

		int[] expectedSizes=new int[] {12947, 12863, 12121, 13008, 13247, 12621, 13515, 13164, 12994, 13728, 12558, 12716, 12922, 13640, 14059, 14519, 13873, 13905, 13586, 14301, 14008, 14000, 13647,
				14254, 14122, 13397, 14887, 13874, 13798, 14489, 14233, 14007, 13566, 15137, 13792, 13668, 14333, 15410, 14122, 14181, 15481, 15281, 14351, 14593, 14310, 14418, 14462, 14785, 14622,
				14687, 14659, 13923, 13929, 13496, 14808, 13492, 13944, 14091, 14558, 14221, 13911, 14424, 15009, 14365, 13643, 14993, 13765, 13756, 14815, 14678, 13506, 14730, 14072, 14664, 14113,
				14842, 12950, 13992, 14175, 14092, 13183, 13616, 13661, 13656, 12513, 13931, 13640, 13776, 14047, 12697, 12775, 13310, 11787, 13537, 13241, 11704, 14093, 12448, 13417, 13056};
		assertEquals(expectedSizes.length, binCounters.length);
		for (int i=0; i<binCounters.length; i++) {
			assertEquals(expectedSizes[i], binCounters[i].size());
		}

		String peptide="ILQEGVDPK";
		double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, (byte)2);
		int index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(20, index);

		double[] keys=binCounters[index].keys();
		Arrays.sort(keys);

		FragmentationModel model=new FragmentationModel(peptide, PARAMETERS.getAAConstants());
		double[] ions=model.getPrimaryIons(PARAMETERS.getFragType());
		int[] expectedCounts=new int[] {181, 388, 34, 17, 10, 16, 10, 14, 10, 10, 5, 14, 5, 11, 6, 5, 15, 18};
		for (int i=0; i<ions.length; i++) {
			double[] matches=PARAMETERS.getFragmentTolerance().getMatches(keys, ions[i]);

			if (matches.length>0) {
				int total=0;
				for (int j=0; j<matches.length; j++) {
					int value=binCounters[index].get(matches[j]);
					total+=value;
				}
				assertEquals(expectedCounts[i], total);
			}
		}

		peptide="FGGGSVELLK";
		mz=PARAMETERS.getAAConstants().getChargedMass(peptide, (byte)2);
		index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(20, index);

		keys=binCounters[index].keys();
		Arrays.sort(keys);

		model=new FragmentationModel(peptide, PARAMETERS.getAAConstants());
		ions=model.getPrimaryIons(PARAMETERS.getFragType());
		expectedCounts=new int[] {388, 41, 8, 60, 6, 3, 9, 3, 10, 14, 7, 8, 6, 7, 10, 6, 15, 22, 29, 32};
		for (int i=0; i<ions.length; i++) {
			double[] matches=PARAMETERS.getFragmentTolerance().getMatches(keys, ions[i]);

			if (matches.length>0) {
				int total=0;
				for (int j=0; j<matches.length; j++) {
					int value=binCounters[index].get(matches[j]);
					total+=value;
				}
				assertEquals(expectedCounts[i], total);
			}
		}

	}
}
