package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Collections;

import com.google.common.collect.Sets;

import junit.framework.TestCase;

/**
 * @author Seth Just
 * @since 12:59 PM 9/22/17
 */
public class FastaPeptideEntryTest extends TestCase {
	/**
	 * Test that new comparison is equivalent to concatenating
	 * accessions in lexicographic order and comparing the
	 * resulting strings.
	 */
	public void testCompareAccessions() throws Exception {
		final FastaPeptideEntry a, b, c, d;

		// shuffled order when creating these to ensure that
		// sorting is correctly applied, if the order were
		// maintained then the assertions below would fail
		// (for example, "2;1" > "1;3;2" but a > b)
		a = of("1", "3", "2"); // 1, 2, 3
		b = of("2", "1");      // 1, 2
		c = of("3", "1");      // 1, 3
		d = of("4", "2", "1"); // 1, 2, 4

		// Test reflexivity
		assertComp(0, a, a);
		assertComp(0, b, b);
		assertComp(0, c, c);
		assertComp(0, d, d);

		// Test lexicographic ordering
		assertComp(1, a, b);
		assertComp(-1, a, c);
		assertComp(-1, a, d);
		assertComp(-1, b, c);
		assertComp(-1, b, d);
		assertComp(1, c, d);

		// Test symmetry
		assertComp(-1, b, a);
		assertComp(1, c, a);
		assertComp(1, d, a);
		assertComp(1, c, b);
		assertComp(1, d, b);
		assertComp(-1, d, c);
	}
	
	public void testModStripping() {
		assertEquals("SSDEENGPPSSPDLDR", new FastaPeptideEntry("S[79.966331]SDEENGPPSSPDLDR").getSequenceWithModsStripped());
		assertEquals("SSDEENGPPSSPDLDR", new FastaPeptideEntry("SS[79.966331]DEENGPPSSPDLDR").getSequenceWithModsStripped());
		assertEquals("SSDEENGPPSSPDLDR", new FastaPeptideEntry("SSDEENGPPS[79.966331]SPDLDR").getSequenceWithModsStripped());
		assertEquals("SSDEENGPPSSPDLDR", new FastaPeptideEntry("SSDEENGPPSS[79.966331]PDLDR").getSequenceWithModsStripped());
		
	}

	private static FastaPeptideEntry of(String... accs) {
		return new FastaPeptideEntry("", Sets.newHashSet(accs), "");
	}

	/**
	 * Copy-paste copy of old method that generated a single
	 * accession from a FastaPeptideEntry.
	 */
	private static String acc(FastaPeptideEntry peptideEntry) {
		ArrayList<String> list = new ArrayList<String>(peptideEntry.getAccessions());
		Collections.sort(list);
		StringBuilder sb = new StringBuilder();
		for (String string : list) {
			if (sb.length() > 0) {
				sb.append(PSMData.ACCESSION_TOKEN);
			}
			sb.append(string);
		}
		return sb.toString();
	}

	private static void assertComp(int expected, FastaPeptideEntry a, FastaPeptideEntry b) {
		assertComp(expected, a.compareTo(b));
		assertComp(a, b);
	}

	private static void assertComp(FastaPeptideEntry a, FastaPeptideEntry b) {
		assertComp(acc(a).compareTo(acc(b)), a.compareTo(b));
	}

	private static void assertComp(int expected, int actual) {
		assertEquals("comparison was wrong; expected: " + expected + ", actual: " + actual, Integer.signum(expected), Integer.signum(actual));
	}
}