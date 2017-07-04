package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import junit.framework.TestCase;

public class FragmentIonBlacklistTest extends TestCase {
	public void testBlacklist() {
		MassTolerance tolerance=new MassTolerance(100);
		FragmentIonBlacklist blacklist=new FragmentIonBlacklist(tolerance);
		
		assertFalse(blacklist.isBlacklisted(100.001, 10f));
		
		blacklist.addIonToBlacklist(100.0, new Range(5f, 15f));
		assertTrue(blacklist.isBlacklisted(100.001, 10f));
		assertFalse(blacklist.isBlacklisted(100.001, 20f));

		blacklist.addIonToBlacklist(100.0, new Range(15f, 25f));
		assertTrue(blacklist.isBlacklisted(100.001, 10f));
		assertTrue(blacklist.isBlacklisted(100.001, 20f));
		assertFalse(blacklist.isBlacklisted(200.001, 20f));
	}

}
