package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
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
		TDoubleIntHashMap[] binCounters=BackgroundGenerator.generateBackground(binArray, entries, true, PARAMETERS).x;

		int[] expectedSizes=new int[] {22304, 22776, 20464, 22246, 22963, 21688, 22642, 22054, 22086, 23288, 21385, 21099, 22227, 22811, 23823, 23851, 22989, 22937, 22540, 24180, 22933, 23294, 22690,
				23721, 23387, 22575, 24529, 22983, 22736, 23230, 23185, 22450, 22756, 24773, 22465, 22563, 23405, 25549, 23162, 22418, 24570, 25482, 22735, 23794, 23093, 23225, 23617, 23360, 24293,
				23101, 23467, 22381, 22140, 21217, 23574, 21759, 22482, 22764, 23304, 23135, 21582, 23212, 23524, 23082, 21996, 23330, 21694, 22175, 23462, 23154, 21503, 23409, 22670, 23195, 23472,
				22812, 20104, 21767, 23572, 21514, 20528, 21396, 22106, 20500, 19482, 23038, 20752, 22490, 21545, 19864, 19827, 20494, 19252, 20311, 20560, 17937, 21410, 19917, 20206, 20786};
		assertEquals(expectedSizes.length, binCounters.length);
		for (int i=0; i<binCounters.length; i++) {
			assertEquals(expectedSizes[i], binCounters[i].size());
		}

		String peptide="ILQEGVDPK";
		byte charge=(byte)2;
		double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
		int index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(19, index);

		double[] keys=binCounters[index].keys();
		Arrays.sort(keys);

		PecanOneFragmentationModel model=new PecanOneFragmentationModel(peptide, PARAMETERS.getAAConstants());
		double[] ions=model.getPrimaryIons(PARAMETERS.getFragType(), charge);
		int[] expectedCounts=new int[] {184, 356, 33, 24, 8, 19, 10, 16, 11, 15, 7, 10, 4, 8, 6, 9, 16, 17};
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
		PecanLibraryEntry entry=model.getPecanSpectrum(charge, keys, binCounters[index], new Range(0f, 200000f), PARAMETERS, false);
		assertEquals(45.135918f, entry.getEuclidianDistance(), 0.0001f);

		peptide="FGGGSVELLK";
		mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
		index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(20, index);

		keys=binCounters[index].keys();
		Arrays.sort(keys);

		model=new PecanOneFragmentationModel(peptide, PARAMETERS.getAAConstants());
		ions=model.getPrimaryIons(PARAMETERS.getFragType(), charge);
		expectedCounts=new int[] {392, 41, 5, 62, 6, 5, 10, 5, 12, 18, 8, 10, 7, 14, 37, 8, 15, 22, 29, 32};
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

		entry=model.getPecanSpectrum(charge, keys, binCounters[index], new Range(0f, 200000f), PARAMETERS, false);
		assertEquals(49.48101f, entry.getEuclidianDistance(), 0.0001f);
	}
}
