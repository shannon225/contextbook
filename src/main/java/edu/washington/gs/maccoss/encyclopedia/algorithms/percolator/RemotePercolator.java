package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;


import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

/**
 * {@link PercolatorVersion} implementation that runs a binary from
 * an arbitrary URI.
 *
 * Subclasses {@link LocalPercolator} to reuse implementations around
 * local paths.
 */
public class RemotePercolator extends LocalPercolator {
	private final URI uri;

	/**
	 * Instantiate a {@code PercolatorVersion} that will run the binary at the specified
	 * URI. Note that the location is not downloaded immediately, but instead only downloaded
	 * when required to run the binary.
	 *
	 * @param uri Not {@code null}
	 * @throws IOException if an I/O error occurs setting up the necessary temporary path
	 */
	public RemotePercolator(URI uri) throws IOException {
		super(Files.createTempFile("Percolator-", ".exe"));
		this.uri = uri;

		final File tmpFile = getFile();

		assert null != tmpFile;
		assert !Files.exists(tmpFile.toPath());

		// Ensure the JVM will clean up the file (eventually).
		// Super method just returns this instance's file.
		tmpFile.deleteOnExit();
	}

	private File getFile() {
		return super.getPercolator();
	}

	private final void downloadPercolator() throws IOException {
		try (InputStream is = uri.toURL().openStream()) {
			Files.copy(is, getFile().toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
	}

	@Override
	public File getPercolator() {
		final File tmpFile = getFile();

		// Check if the file's already been downloaded.
		if (!Files.exists(tmpFile.toPath())) {
			try {
				downloadPercolator();
			} catch (IOException e) {
				throw new UncheckedIOException(e);
			}
		}

		return tmpFile;
	}

	@Override
	public String toString() {
		return uri.toString();
	}
}
