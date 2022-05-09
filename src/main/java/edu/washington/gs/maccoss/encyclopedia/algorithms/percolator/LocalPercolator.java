package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Implementation for a Percolator version that exists as a local file.
 */
public class LocalPercolator extends ExternalPercolator {
	/**
	 * Nonnull
	 */
	private final Path exe;

	/**
	 * @param exe Not {@code null}
	 */
	public LocalPercolator(File exe) {
		this(exe.toPath());
	}

	/**
	 * @param exe Not {@code null}
	 */
	public LocalPercolator(Path exe) {
		this.exe = exe;

		if (!Files.isExecutable(exe)) {
			Logger.errorLine("Percolator file is not executable; continuing anyway, though this may cause errors. Location: " + exe.toAbsolutePath());
		}
	}

	@Override
	public File getPercolator() {
		return exe.toFile();
	}
}
