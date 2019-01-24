package edu.washington.gs.maccoss.encyclopedia.utils;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import junit.framework.TestCase;

public class VersioningDetectorTest extends TestCase {
	public void testVersionDetection() {
		boolean pass=VersioningDetector.checkVersionGUI(ProgramType.EncyclopeDIA, null, true, false);
		assertFalse(pass);
	}

}
