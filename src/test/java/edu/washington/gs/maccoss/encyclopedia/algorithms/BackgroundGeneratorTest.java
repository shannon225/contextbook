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
	private static final SearchParameters PARAMETERS=new SearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50), DigestionEnzyme.getEnzyme("trypsin"));

	public void testGenerateBackground() {
		InputStream is=getClass().getResourceAsStream("/mouse_20150911_uniprot_sp.fasta");
		ArrayList<FastaEntry> entries=FastaReader.readFasta(is, "mouse_20150911_uniprot_sp.fasta");
		TDoubleArrayList bins=new TDoubleArrayList();
		for (double i=400.0; i<=900.0; i+=5.0) {
			bins.add(i);
		}
		double[] binArray=bins.toArray();
		TDoubleIntHashMap[] binCounters=BackgroundGenerator.generateBackground(binArray, entries, PARAMETERS).x;

		int[] expectedSizes=new int[] {61931, 63402, 64018, 65614, 66799, 67520, 68660, 69868, 70879, 71698, 72725, 74656, 74381, 75946, 75106, 76099, 77119, 78127, 78183, 79101, 78580, 80176, 81050,
				81264, 81261, 82083, 81800, 82215, 83790, 83489, 84408, 84262, 84977, 86106, 85042, 86042, 84445, 85417, 84377, 86291, 84028, 84687, 87367, 85448, 83981, 86387, 86753, 85877, 86506,
				87500, 85405, 87050, 86628, 88037, 83770, 86467, 86291, 88143, 85045, 87602, 87065, 87145, 88023, 87734, 88110, 88349, 86228, 84490, 84683, 83114, 83774, 85508, 85082, 83354, 85966,
				83674, 83093, 82152, 85346, 81338, 82518, 80409, 79875, 86097, 81132, 81157, 80317, 82413, 81383, 81660, 83528, 78853, 81202, 79133, 78277, 82234, 82806, 77735, 78448, 81221};
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
		int[] expectedCounts=new int[] {1771, 5280, 353, 318, 97, 231, 90, 125, 93, 74, 79, 109, 162, 123, 189, 136, 322, 320};
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
		expectedCounts=new int[] {1771, 5280, 353, 318, 97, 231, 90, 125, 93, 74, 79, 109, 162, 123, 189, 136, 322, 320};
		for (int i=0; i<ions.length; i++) {
			double[] matches=PARAMETERS.getFragmentTolerance().getMatches(keys, ions[i]);

			if (matches.length>0) {
				int total=0;
				for (int j=0; j<matches.length; j++) {
					int value=binCounters[index].get(matches[j]);
					total+=value;
				}
				System.out.println(ions[i]+"\t"+total);
				// assertEquals(expectedCounts[i], total);
			}
		}

	}
}
