package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.HashMap;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.ModificationMassMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import gnu.trove.map.hash.TCharDoubleHashMap;
import junit.framework.TestCase;

public class DotProductTest extends TestCase {
	
	public void testEmptyDotProduct() {
		LibraryEntry entry=getEntry(new double[] {}, new float[] {});
		FragmentScan spectrum=getStripe(new double[] {}, new float[] {});
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-ftol", "100");
		map.put("-lftol", "100");
		assertEquals(0.0f, new DotProduct(SearchParameterParser.parseParameters(map)).score(entry, spectrum, new float[] {}, null));
	}
	public void testDotProduct() {
		LibraryEntry entry=getEntry(new double[] {1.0, 29.0, 300.01, 1000.0, 1200.0}, new float[] {7, 7, 2, 3, 7});
		FragmentScan spectrum=getStripe(new double[] {30.0, 300.0, 1001.0, 1300.0}, new float[] {7, 10, 4, 7});
		
		HashMap<String, String> map=SearchParameterParser.getDefaultParameters();
		map.put("-ftol", "1000");
		map.put("-lftol", "1000");
		assertEquals(32.0f, new DotProduct(SearchParameterParser.parseParameters(map)).score(entry, spectrum, new float[] {}, null));

		map.put("-ftol", "100");
		map.put("-lftol", "100");
		assertEquals(20.0f, new DotProduct(SearchParameterParser.parseParameters(map)).score(entry, spectrum, new float[] {}, null));

		map.put("-ftol", "10");
		map.put("-lftol", "10");
		assertEquals(0.0f, new DotProduct(SearchParameterParser.parseParameters(map)).score(entry, spectrum, new float[] {}, null));
	}
	
	public LibraryEntry getEntry(double[] masses, float[] intensities) {
		return new LibraryEntry("", new HashSet<String>(), 1, (byte)1, "", 1, 1, 1, masses, intensities, new AminoAcidConstants(new TCharDoubleHashMap(), new ModificationMassMap()));
	}
	public FragmentScan getStripe(double[] masses, float[] intensities) {
		return new FragmentScan("", "", 1, 1, 0, 0f, 1, 1, masses, intensities);
	}
}
