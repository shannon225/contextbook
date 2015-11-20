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
		InputStream is=getClass().getResourceAsStream("/mouse_20150911_uniprot_sp.fasta");
		ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "mouse_20150911_uniprot_sp.fasta");
		TDoubleArrayList bins=new TDoubleArrayList();
		for (double i=400.0; i<=900.0; i+=5.0) {
			bins.add(i);
		}
		double[] binArray=bins.toArray();
		TDoubleIntHashMap[] binCounters=BackgroundGenerator.generateBackground(binArray, entries, PARAMETERS).x;

		int[] expectedSizes=new int[] {62490, 63036, 64157, 66232, 66727, 66942, 68383, 68992, 70288, 71440, 72849, 73478, 74693, 76248, 74334, 76749, 77078, 78343, 78905, 78839, 78156, 80610, 80232,
				80911, 80508, 82407, 81740, 81545, 83369, 82531, 83396, 82540, 84143, 85849, 85593, 85497, 85755, 85739, 85708, 86573, 84195, 84338, 85280, 86974, 84606, 86134, 86992, 86090, 86174,
				86553, 84988, 85015, 85303, 88357, 84788, 87080, 86404, 87924, 85198, 86969, 85578, 86588, 86919, 88670, 88348, 88526, 85224, 85468, 85437, 83064, 83822, 85938, 85993, 82934, 85597,
				84000, 81408, 81229, 84133, 81342, 82883, 79410, 80171, 84549, 78541, 81222, 81128, 80980, 80707, 82055, 82072, 79129, 81366, 78730, 79943, 81257, 80397, 78726, 79911, 80636};
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
		int[] expectedCounts=new int[] {1704, 5252, 323, 328, 85, 232, 73, 134, 84, 73, 85, 98, 156, 112, 185, 126, 304, 300};
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
		expectedCounts=new int[] {5252, 463, 44, 903, 52, 8, 158, 43, 69, 122, 76, 63, 43, 100, 144, 94, 131, 244, 300, 309};
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
