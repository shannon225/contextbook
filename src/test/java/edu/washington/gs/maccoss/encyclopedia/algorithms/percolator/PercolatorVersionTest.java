package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;
import org.junit.AssumptionViolatedException;
import org.junit.Ignore;
import org.junit.Test;

import java.net.URI;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class PercolatorVersionTest {
	@Test
	public void testParsePercolator() throws Exception {
		assertEquals(PercolatorVersion.v2p10, PercolatorVersion.getVersion(PercolatorVersion.V2_10));
		assertEquals(PercolatorVersion.v3p01, PercolatorVersion.getVersion(PercolatorVersion.V3_01));
		assertEquals(PercolatorVersion.v3p05, PercolatorVersion.getVersion(PercolatorVersion.V3_05));
		assertEquals(PercolatorVersion.v3p01, PercolatorVersion.getVersion("3"));
		assertEquals(PercolatorVersion.v3p05, PercolatorVersion.getVersion("3.5"));
		assertEquals(PercolatorVersion.v3p01, PercolatorVersion.getVersion("3.1"));

		// local path
		PercolatorVersion parsed = PercolatorVersion.getVersion(PercolatorVersion.v3p05.getPercolator().getAbsolutePath());
		assertTrue("Did not get expected PercolatorVersion impl!", parsed instanceof LocalPercolator);
		assertEquals(PercolatorVersion.v3p05.getMajorVersion(), parsed.getMajorVersion());
	}

	@Test
	public void testParsePercolatorLocalURI() throws Exception {
		PercolatorVersion parsed = PercolatorVersion.getVersion(PercolatorVersion.v2p10.getPercolator().toURI().toString());
		assertTrue("Did not get expected PercolatorVersion impl!", parsed instanceof LocalPercolator);
		assertTrue("Got unexpected PercolatorVersion impl for local URI!", !(parsed instanceof RemotePercolator));
		assertEquals(PercolatorVersion.v2p10.getMajorVersion(), parsed.getMajorVersion());
	}

	@Test
	public void testParsePercolatorURI() throws Exception {
		String uri;
		switch (OSDetector.getOS()) {
			case WINDOWS:
				uri = "https://bitbucket.org/searleb/encyclopedia/raw/encyclopedia-1.12.31/src/main/resources/bin/percolator-v3-01.exe";
				break;
			case MAC:
				uri = "https://bitbucket.org/searleb/encyclopedia/raw/encyclopedia-1.12.31/src/main/resources/bin/percolator-v3-01.mac";
				break;
			case LINUX:
				uri = "https://bitbucket.org/searleb/encyclopedia/raw/encyclopedia-1.12.31/src/main/resources/bin/percolator-v3-01.lin";
				break;
			default:
				throw new AssumptionViolatedException("Can't run test without recognizing the OS!");
		}

		PercolatorVersion parsed = PercolatorVersion.getVersion(uri);
		assertTrue("Did not get expected PercolatorVersion impl!", parsed instanceof ExternalPercolator);
		assertEquals(3, parsed.getMajorVersion());
		assertEquals("3.01", PercolatorExecutor.checkPercolatorVersion(parsed));
	}
}