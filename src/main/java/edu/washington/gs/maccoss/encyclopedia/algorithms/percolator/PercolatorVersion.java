package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.OSDetector;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

public interface PercolatorVersion {
	PercolatorVersion v2p10 = InternalPercolatorVersion.v2p10;
	PercolatorVersion v3p01 = InternalPercolatorVersion.v3p01;
	PercolatorVersion v3p05 = InternalPercolatorVersion.v3p05;

	/**
	 * By default, use version 3. Until we default to Percolator 4,
	 * don't change this definition and instead update {@link #DEFAULT_VERSION_3}.
	 */
	PercolatorVersion DEFAULT_VERSION = PercolatorVersion.DEFAULT_VERSION_3;

	/**
	 * The version of Percolator 3 that should be used by default,
	 * e.g. when running with `-percolatorVersion 3`. Currently 3.01
	 * due to issues observed with 3.05.
	 */
	PercolatorVersion DEFAULT_VERSION_3 = PercolatorVersion.v3p01;

	PercolatorVersion[] VALID_VERSIONS = new PercolatorVersion[]{v3p01, v2p10};

	String V3_05 = "v3-05";
	String V3_01 = "v3-01";
	String V2_10 = "v2-10";

	int UNKNOWN_MAJOR_VERSION = -1;

	static PercolatorVersion getVersion(String s) {
		if (s == null || s.length() == 0) return DEFAULT_VERSION;
		if (V2_10.equals(s)) return v2p10;
		if (V3_01.equals(s)) return v3p01;
		if (V3_05.equals(s)) return v3p05;
		if ("2".equals(s)) return v2p10;
		if ("3".equals(s)) return DEFAULT_VERSION_3;
		if ("2.10".equals(s)) return v2p10;
		if ("3.1".equals(s)) return v3p01;
		if ("3.5".equals(s)) return v3p05;
		if ("3.01".equals(s)) return v3p01;
		if ("3.05".equals(s)) return v3p05;

		final PercolatorVersion parsedAsFile = parseFilePath(s);
		if (null != parsedAsFile) {
			return parsedAsFile;
		}

		final PercolatorVersion parsedAsUri = parseUri(s);
		if (null != parsedAsUri) {
			return parsedAsUri;
		}

		Logger.errorLine("Could not parse Percolator version \"" + s + "\"; falling back to " + DEFAULT_VERSION);

		return DEFAULT_VERSION;
	}

	/**
	 * Check if the given string is a valid filepath, and if so,
	 * return a {@link LocalPercolator} instance that will use it.
	 *
	 * @return a suitable {@code PercolatorVersion}, or {@code null}
	 */
	static PercolatorVersion parseFilePath(String s) {
		if (null == s || s.isEmpty()) {
			return null;
		}

		final Path parsed = Paths.get(s);
		if (Files.exists(parsed)) {
			return new LocalPercolator(parsed);
		}

		return null;
	}

	/**
	 * Check if the given string is a valid URI, and if so,
	 * return a suitable {@code PercolatorVersion} instance that will use it.
	 *
	 * @return a suitable {@code PercolatorVersion}, or {@code null}
	 */
	static PercolatorVersion parseUri(String s) {
		if (null == s || s.isEmpty()) {
			return null;
		}

		final URI parsed;
		try {
			parsed = new URI(s);
		} catch (URISyntaxException e) {
			// Don't bother logging here; the caller will log if the entire
			// parsing process fails.
			return null;
		}

		if ("file".equals(parsed.getScheme())) {
			return new LocalPercolator(Paths.get(parsed.getPath()));
		}

		return null;
	}

	int getMajorVersion();

	/**
	 * Get the local executable file for this version of Percolator, possibly after copying it to temp.
	 */
	File getPercolator();

	enum InternalPercolatorVersion implements PercolatorVersion {
		v2p10, v3p01, v3p05;

		@Override
		public String toString() {
			switch (this) {
				case v2p10:
					return V2_10;
				case v3p01:
					return V3_01;
				case v3p05:
					return V3_05;
				default:
					return DEFAULT_VERSION.toString();
			}
		}

		@Override
		public int getMajorVersion() {
			switch (this) {
				case v2p10:
					return 2;
				case v3p01:
				case v3p05:
				default:
					return 3;
			}
		}

		@Override
		public File getPercolator() {
			try {
				File percolator = File.createTempFile("Percolator-" + this + "-", ".exe");
				percolator.deleteOnExit();

				OSDetector.OS os = OSDetector.getOS();
				switch (os) {
					case WINDOWS: {
						InputStream is = PercolatorExecutor.class.getResourceAsStream("/bin/percolator-" + this + ".exe");
						Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
						percolator.setExecutable(true);

						// not necessary for the crux version of percolator
						//loadLibraryFile(percolator, "xerces-c_3_1.dll");
						//loadLibraryFile(percolator, "msvcr120.dll");
						//loadLibraryFile(percolator, "msvcp120.dll");

						return percolator;
					}
					case MAC: {
						InputStream is = PercolatorExecutor.class.getResourceAsStream("/bin/percolator-" + this + ".mac");
						Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
						percolator.setExecutable(true);
						return percolator;
					}
					case LINUX:
						InputStream is = PercolatorExecutor.class.getResourceAsStream("/bin/percolator-" + this + ".lin");
						Files.copy(is, percolator.toPath(), StandardCopyOption.REPLACE_EXISTING);
						percolator.setExecutable(true);
						return percolator;
				}
				throw new EncyclopediaException("Sorry, Percolator for " + OSDetector.getOSName(os) + " is not set up yet!");
			} catch (IOException ioe) {
				throw new EncyclopediaException("Unexpected exception finding Percolator", ioe);
			}
		}

		static void loadLibraryFile(File percolator, String target) throws IOException {
			File file=new File(percolator.getParentFile(), target);
			file.deleteOnExit();
			InputStream is=PercolatorExecutor.class.getResourceAsStream("/bin/"+target);
			Files.copy(is, file.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
	}
}