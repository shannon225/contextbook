package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanOneFragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaEntryInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import junit.framework.TestCase;

public class BackgroundGeneratorTest extends TestCase {
	private static final PecanSearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"), 0);

	public void testGenerateBackground() {
		InputStream is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(is, "contam_correctNL.fasta");
		TDoubleArrayList bins=new TDoubleArrayList();
		for (double i=400.0; i<=900.0; i+=5.0) {
			bins.add(i);
		}
		double[] binArray=bins.toArray();
		TDoubleIntHashMap[] binCounters=BackgroundGenerator.generateBackground(binArray, entries, true, PARAMETERS).x;

		int[] expectedSizes=new int[] { 22312, 22778, 20460, 22258, 22958, 21681, 22643, 22048, 22089, 23289, 21389, 21103, 22233, 22812, 23816, 23849, 22989, 22942, 22543, 24179, 22925, 23298, 22690,
				23727, 23387, 22575, 24528, 22986, 22739, 23232, 23179, 22451, 22762, 24768, 22461, 22562, 23410, 25550, 23169, 22415, 24572, 25486, 22734, 23794, 23094, 23222, 23623, 23369, 24293,
				23101, 23468, 22380, 22140, 21215, 23576, 21763, 22484, 22764, 23308, 23133, 21584, 23214, 23521, 23085, 21994, 23329, 21697, 22176, 23465, 23156, 21497, 23409, 22671, 23197, 23471,
				22813, 20105, 21768, 23571, 21514, 20512, 21399, 22106, 20500, 19480, 23036, 20754, 22494, 21547, 19864, 19827, 20495, 19250, 20312, 20563, 17937, 21411, 19921, 20206, 20786 };
		assertEquals(expectedSizes.length, binCounters.length);
		for (int i=0; i<binCounters.length; i++) {
			assertEquals(expectedSizes[i], binCounters[i].size());
		}

		String peptide="ILQEGVDPK";
		byte charge=(byte) 2;
		double mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
		int index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(19, index);

		double[] keys=binCounters[index].keys();
		Arrays.sort(keys);

		PecanOneFragmentationModel model=new PecanOneFragmentationModel(new FastaPeptideEntry(peptide), PARAMETERS.getAAConstants());
		double[] ions=model.getPrimaryIons(PARAMETERS.getFragType(), charge);
		int[] expectedCounts=new int[] { 184, 356, 33, 24, 8, 19, 10, 16, 11, 15, 7, 10, 4, 8, 6, 9, 16, 17 };
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
		assertEquals(39.129665f, entry.getEuclidianDistance(), 0.0001f);

		peptide="FGGGSVELLK";
		mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
		index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(20, index);

		keys=binCounters[index].keys();
		Arrays.sort(keys);

		model=new PecanOneFragmentationModel(new FastaPeptideEntry(peptide), PARAMETERS.getAAConstants());
		ions=model.getPrimaryIons(PARAMETERS.getFragType(), charge);
		expectedCounts=new int[] { 392, 41, 5, 62, 6, 5, 10, 5, 12, 18, 8, 10, 7, 14, 37, 8, 15, 22, 29, 32 };
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
		assertEquals(42.785778f, entry.getEuclidianDistance(), 0.0001f);
	}
}
