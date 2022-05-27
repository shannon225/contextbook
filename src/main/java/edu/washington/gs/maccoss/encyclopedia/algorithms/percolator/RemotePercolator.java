package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;


import com.google.common.collect.ImmutableSet;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.PosixFilePermission;

/**
 * {@link PercolatorVersion} implementation that runs a binary from
 * an arbitrary URI.
 *
 * Subclasses {@link LocalPercolator} to reuse implementations around
 * local paths.
 */
public class RemotePercolator extends LocalPercolator {
	private final URI uri;

	private boolean didDownload = false;

	/**
	 * Instantiate a {@code PercolatorVersion} that will run the binary at the specified
	 * URI. Note that the location is not downloaded immediately, but instead only downloaded
	 * when required to run the binary.
	 *
	 * @param uri Not {@code null}
	 * @throws IOException if an I/O error occurs setting up the necessary temporary path
	 */
	public RemotePercolator(URI uri) throws IOException {
		super(setupTempFile());
		this.uri = uri;
	}

	private static Path setupTempFile() throws IOException {
		final Path tmpFile = Files.createTempFile("Percolator-", ".exe");

		// Ensure the JVM will clean up the file (eventually).
		// Super method just returns this instance's file.
		tmpFile.toFile().deleteOnExit();

		final boolean setEx = tmpFile.toFile().setExecutable(true);

		assert setEx; // This line runs only during tests, or if java's invoked with -ea

		return tmpFile;
	}

	private File getFile() {
		return super.getPercolator();
	}

	private void downloadPercolator() throws IOException {
		final Path path = getFile().toPath();

		try (InputStream is = uri.toURL().openStream()) {
			Files.copy(is, path, StandardCopyOption.REPLACE_EXISTING);
		}

		if (!path.toFile().setExecutable(true)) {
			throw new IOException("Could not set executable flag on downloaded Percolator!");
		}
	}

	@Override
	public File getPercolator() {
		final File tmpFile = getFile();

		// Check if the file's already been downloaded.
		if (!didDownload) {
			try {
				downloadPercolator();
			} catch (IOException e) {
				throw new UncheckedIOException(e);
			}
			didDownload = true; // only set if we succeeded; in rare cases this could lead to repeated retries!
		}

		return tmpFile;
	}

	@Override
	public String toString() {
		return uri.toString();
	}
}
