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
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.FastaReader;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TDoubleIntHashMap;
import junit.framework.TestCase;

public class BackgroundGeneratorTest extends TestCase {
	private static final PecanSearchParameters PARAMETERS=new PecanSearchParameters(new AminoAcidConstants(), FragmentationType.CID, new MassTolerance(50), new MassTolerance(50),
			DigestionEnzyme.getEnzyme("trypsin"), 0, false, true, false);

	public void testGenerateBackground() {
		InputStream is=getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta");
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		ArrayList<FastaEntryInterface> entries=FastaReader.readFasta(is, "contam_correctNL.fasta", parameters);
		TDoubleArrayList bins=new TDoubleArrayList();
		for (double i=400.0; i<=900.0; i+=5.0) {
			bins.add(i);
		}
		double[] binArray=bins.toArray();
		TDoubleIntHashMap[] binCounters=BackgroundGenerator.generateBackground(binArray, entries, true, PARAMETERS).x;

		int[] expectedSizes=new int[] { 22641, 23383, 20904, 22813, 23452, 22317, 23223, 22302, 22625, 23667, 21814,
				21687, 22903, 23534, 24418, 24216, 23604, 23444, 23108, 24814, 23620, 23785, 23464, 24267, 23922, 22855,
				25125, 23665, 23179, 23557, 23622, 23033, 23077, 25219, 22846, 23151, 24196, 25932, 23860, 23205, 24936,
				25945, 23355, 24353, 23507, 23673, 24012, 23788, 24843, 23623, 24000, 22809, 22662, 21927, 24152, 22270,
				23175, 23293, 23851, 23871, 22192, 23572, 24249, 23757, 22757, 23954, 21994, 22740, 23775, 23962, 22179,
				23659, 23327, 23508, 24177, 23386, 21067, 22488, 24014, 22120, 21101, 22136, 22816, 21035, 19846, 23885,
				21121, 22843, 21679, 20857, 20605, 21223, 19493, 20922, 21385, 18352, 21979, 20488, 20583, 21308
		};
		assertEquals(expectedSizes.length, binCounters.length);
		for (int i=0; i<binCounters.length; i++) {
			assertEquals(expectedSizes[i], binCounters[i].size(), 200);
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
		double[] ions=model.getPrimaryIons(PARAMETERS.getFragType(), charge, false);
		int[] expectedCounts=new int[] { 190, 368, 36, 25, 8, 20, 11, 17, 11, 15, 7, 10, 4, 8, 6, 9, 17, 18 };
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
		assertEquals(38.76069f, entry.getEuclidianDistance(), 0.0001f);

		peptide="FGGGSVELLK";
		mz=PARAMETERS.getAAConstants().getChargedMass(peptide, charge);
		index=Arrays.binarySearch(binArray, mz);
		index=(-(index+1))-1;
		assertEquals(20, index);

		keys=binCounters[index].keys();
		Arrays.sort(keys);

		model=new PecanOneFragmentationModel(new FastaPeptideEntry(peptide), PARAMETERS.getAAConstants());
		ions=model.getPrimaryIons(PARAMETERS.getFragType(), charge, false);
		expectedCounts=new int[] { 406, 43, 5, 65, 6, 5, 11, 5, 12, 18, 8, 11, 7, 14, 38, 8, 15, 23, 29, 32 };
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
		assertEquals(42.445232f, entry.getEuclidianDistance(), 0.0001f);
	}
}
