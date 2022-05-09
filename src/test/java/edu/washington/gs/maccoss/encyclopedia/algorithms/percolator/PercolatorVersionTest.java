package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import org.junit.Test;

import static org.junit.Assert.*;

public class PercolatorVersionTest {
	@Test
	public void testParsePercolator() throws Exception {
		assertEquals(PercolatorVersion.v2p10, PercolatorVersion.getVersion(PercolatorVersion.V2_10));
		assertEquals(PercolatorVersion.v3p01, PercolatorVersion.getVersion(PercolatorVersion.V3_01));
		assertEquals(PercolatorVersion.v3p05, PercolatorVersion.getVersion(PercolatorVersion.V3_05));
		assertEquals(PercolatorVersion.v3p05, PercolatorVersion.getVersion("3"));
		assertEquals(PercolatorVersion.v3p05, PercolatorVersion.getVersion("3.5"));
		assertEquals(PercolatorVersion.v3p01, PercolatorVersion.getVersion("3.1"));

		// Parse local path
		PercolatorVersion parsed = PercolatorVersion.getVersion(PercolatorVersion.v3p05.getPercolator().getAbsolutePath());
		assertTrue(parsed instanceof ExternalPercolator);
		assertEquals(PercolatorVersion.v3p05.getMajorVersion(), parsed.getMajorVersion());
	}
}