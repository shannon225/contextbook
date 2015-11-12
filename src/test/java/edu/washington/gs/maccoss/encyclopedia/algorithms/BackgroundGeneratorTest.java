package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import junit.framework.TestCase;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaEntry;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;

public class BackgroundGeneratorTest extends TestCase {
	private static final SearchParameters PARAMETERS=new SearchParameters(FragmentationType.CID, new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"));

	public void testGenerateBackground() {
		InputStream is=getClass().getResourceAsStream("/mouse_UP000000589_10090.fasta");
		ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "mouse_UP000000589_10090.fasta");
		TDoubleArrayList bins=new TDoubleArrayList();
		for (double i=400.0; i<=900.0; i+=5.0) {
			bins.add(i);
		}
		double[] binArray=bins.toArray();
		TDoubleIntHashMap[] binCounters=BackgroundGenerator.generateBackground(binArray, entries, PARAMETERS);

		int[] expectedSizes=new int[] {68095, 70163, 71351, 72651, 74186, 74848, 76046, 78188, 78463, 80511, 81045, 83142, 82816, 84382, 84005, 85045, 86437, 87761, 87677, 88600, 88840, 89836, 90971,
				91858, 91891, 92636, 92624, 92373, 93950, 93511, 96287, 95370, 96097, 97730, 96298, 96954, 95261, 95809, 95025, 96885, 95036, 95583, 97803, 96567, 95511, 97193, 97607, 97554, 96753,
				98948, 97446, 99369, 98232, 100442, 95683, 98914, 97806, 99973, 97097, 98354, 98220, 98915, 99487, 99428, 98645, 99585, 97484, 95589, 95668, 94619, 96268, 97633, 96619, 95450, 97024,
				96701, 95202, 94360, 98308, 93997, 96584, 92989, 91915, 99797, 93586, 93947, 93504, 94312, 93228, 93236, 95676, 90824, 93096, 91591, 89875, 93242, 94149, 88741, 90810, 92255};
		assertEquals(expectedSizes.length, binCounters.length);
		for (int i=0; i<binCounters.length; i++) {
			assertEquals(expectedSizes[i], binCounters[i].size());
		}

		String peptide="ILQEGVDPK";
		double mz=MassConstants.getChargedMass(peptide, (byte)2);
		int index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(20, index);

		double[] keys=binCounters[index].keys();
		Arrays.sort(keys);

		FragmentationModel model=new FragmentationModel(peptide);
		double[] ions=model.getPrimaryIons(PARAMETERS.getFragType());
		int[] expectedCounts=new int[] {2007, 5572, 377, 356, 124, 178, 72, 56, 92, 115, 30, 40, 26, 27, 14, 17, 6, 16};
		for (int i=0; i<ions.length; i++) {
			double[] matches=PARAMETERS.getTolerance().getMatches(keys, ions[i]);

			if (matches.length>0) {
				int total=0;
				for (int j=0; j<matches.length; j++) {
					total+=binCounters[0].get(matches[j]);
				}
				assertEquals(expectedCounts[i], total);
			}
		}
	}
}
