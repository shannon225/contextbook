package edu.washington.gs.maccoss.encyclopedia.filereaders;

import junit.framework.TestCase;

public class SpectronautCSVToLibraryConverterTest extends TestCase {
	public void testParseMods() {
		String seq=SpectronautCSVToLibraryConverter.parseMods("_YNDLEKHVCSIK_");
		assertEquals("YNDLEKHVCSIK", seq);
		seq=SpectronautCSVToLibraryConverter.parseMods("_YNDLEKHVC[Carbamidomethyl (C)]SIK_");
		assertEquals("YNDLEKHVC[+57.021464]SIK", seq);
		seq=SpectronautCSVToLibraryConverter.parseMods("_[Acetyl (Protein N-term)]YNDLEKHVCSIK_");
		assertEquals("Y[+42.010565]NDLEKHVCSIK", seq);
	}
}
