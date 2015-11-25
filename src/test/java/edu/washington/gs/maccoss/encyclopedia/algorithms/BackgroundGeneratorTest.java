package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import junit.framework.TestCase;

public class BackgroundGeneratorTest extends TestCase {
	private static final SearchParameters PARAMETERS=new SearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"), 0);

	public void testGenerateBackground() {
		InputStream is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "contam_correctNL.fasta");
		TDoubleArrayList bins=new TDoubleArrayList();
		for (double i=400.0; i<=900.0; i+=5.0) {
			bins.add(i);
		}
		double[] binArray=bins.toArray();
		TDoubleIntHashMap[] binCounters=BackgroundGenerator.generateBackground(binArray, entries, PARAMETERS).x;

		int[] expectedSizes=new int[] {22411, 22285, 20803, 22091, 22989, 21518, 22927, 21994, 21927, 23646, 21077, 21235, 21896, 22808, 23758, 24460, 22840, 22963, 22534, 24042, 22939, 23248, 22623,
				23713, 23381, 22541, 24799, 22672, 22485, 23836, 23091, 22896, 22330, 24781, 22344, 22326, 23640, 25040, 23303, 22636, 24878, 24936, 22913, 23791, 23293, 23190, 22971, 23927, 23546,
				23751, 23735, 22263, 22247, 21292, 23214, 21628, 22859, 22587, 23317, 22818, 21889, 22941, 23696, 22757, 21865, 23807, 22041, 21942, 23367, 23338, 21370, 23542, 22545, 23379, 22672,
				23597, 20143, 21939, 22773, 22119, 20293, 21422, 22269, 20827, 19504, 22253, 21577, 22136, 21776, 19706, 19647, 20949, 18634, 21040, 20212, 17847, 21665, 19123, 20758, 21061};
		assertEquals(expectedSizes.length, binCounters.length);
		for (int i=0; i<binCounters.length; i++) {
			assertEquals(expectedSizes[i], binCounters[i].size());
		}

		String peptide="ILQEGVDPK";
		byte charge=(byte)2;
		double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
		int index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(20, index);

		double[] keys=binCounters[index].keys();
		Arrays.sort(keys);

		PecanOneFragmentationModel model=new PecanOneFragmentationModel(peptide, PARAMETERS.getAAConstants());
		double[] ions=model.getPrimaryIons(PARAMETERS.getFragType(), charge);
		int[] expectedCounts=new int[] {198, 388, 34, 20, 12, 21, 12, 16, 12, 16, 6, 16, 25, 11, 6, 5, 15, 18};
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
		PecanLibraryEntry entry=model.getPecanSpectrum(charge, keys, binCounters[index], PARAMETERS, false);
		assertEquals(34.20645f, entry.getEuclidianDistance(), 0.0001f);

		peptide="FGGGSVELLK";
		mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
		index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(20, index);

		keys=binCounters[index].keys();
		Arrays.sort(keys);

		model=new PecanOneFragmentationModel(peptide, PARAMETERS.getAAConstants());
		ions=model.getPrimaryIons(PARAMETERS.getFragType(), charge);
		expectedCounts=new int[] {388, 41, 8, 60, 7, 5, 9, 5, 12, 17, 7, 8, 6, 14, 33, 6, 15, 22, 29, 32};
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

		entry=model.getPecanSpectrum(charge, keys, binCounters[index], PARAMETERS, false);
		assertEquals(43.045853f, entry.getEuclidianDistance(), 0.0001f);
	}
}
