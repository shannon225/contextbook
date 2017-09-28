package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import junit.framework.TestCase;

import java.io.File;

public class EncyclopediaJobDataTest extends TestCase {
	public void testGetOutputAbsolutePathPrefix() {
		String absolutePath="/Users/searleb/Documents/projects/encyclopedia/encyclopedia/121115_BCS_HeLa_24mz_400_1000.dia.encyclopedia.txt";
		assertEquals("/Users/searleb/Documents/projects/encyclopedia/encyclopedia/121115_BCS_HeLa_24mz_400_1000.dia", EncyclopediaJobData.getPrefixFromOutput(new File(absolutePath)));

		absolutePath="/Users/searleb/Documents/projects/encyclopedia/encyclopedia/121115_BCS_HeLa_24mz_400_1000.dia.encyclopedia";
		assertEquals("/Users/searleb/Documents/projects/encyclopedia/encyclopedia/121115_BCS_HeLa_24mz_400_1000.dia.encyclopedia", EncyclopediaJobData.getPrefixFromOutput(new File(absolutePath)));
	}

}
