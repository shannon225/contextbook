package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import org.junit.Test;

import java.io.File;

import static org.junit.Assert.*;

public class ExternalPercolatorTest {

	@Test
	public void getMajorVersion() {
		checkMajorVersion(PercolatorVersion.v3p05);
		checkMajorVersion(PercolatorVersion.v3p01);
		checkMajorVersion(PercolatorVersion.v2p10);
	}

	static void checkMajorVersion(PercolatorVersion version) {
		assertEquals(version.getMajorVersion(), externalGetMajorVersion(version));
	}

	/**
	 * Wrap the given version in an {@link ExternalPercolator} and
	 * return its parsed major version.
	 */
	static int externalGetMajorVersion(PercolatorVersion version) {
		final ExternalPercolator perc = new ExternalPercolator() {
			@Override
			public File getPercolator() {
				return version.getPercolator();
			}
		};

		return perc.getMajorVersion();
	}
}