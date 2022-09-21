package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Base class for non-builtin Percolator versions.
 */
public abstract class ExternalPercolator implements PercolatorVersion {
	// anchored only at start, so suffixes are ignored
	private static final Pattern VERSION_PATTERN = Pattern.compile("^(\\d+)(?:\\.(\\d+)(?:\\.(\\d+))?)?");

	@Override
	public int getMajorVersion() {
		final String version;
		try {
			version = PercolatorExecutor.checkPercolatorVersion(this);
		} catch (IOException | InterruptedException e) {
			int defaultVersion = PercolatorVersion.DEFAULT_VERSION.getMajorVersion();
			Logger.errorLine("Unable to check Percolator version; assuming default " + defaultVersion);
			return defaultVersion;
		}

		final Matcher m = VERSION_PATTERN.matcher(version);

		if (m.find()) {
			return Integer.parseInt(m.group(1));
		}
		return UNKNOWN_MAJOR_VERSION;
	}
}
