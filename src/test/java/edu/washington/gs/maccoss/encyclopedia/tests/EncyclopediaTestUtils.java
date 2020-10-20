package edu.washington.gs.maccoss.encyclopedia.tests;

import org.apache.commons.io.FileUtils;
import org.junit.Assume;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

public abstract class EncyclopediaTestUtils {
	private EncyclopediaTestUtils() throws IllegalAccessError {
		throw new IllegalAccessError();
	}

	public static Path getResourceAsTempFile(Class<?> klass, String resource, Path tmpDir, String name, String suffix) throws IOException {
		final Path dest = Files.createTempFile(tmpDir, name, suffix);
		FileUtils.forceDeleteOnExit(dest.toFile());
		copyResourceToFile(klass, resource, dest);
		return dest;
	}

	public static void copyResourceToFile(Class<?> klass, String resource, Path dest) throws IOException {
		try (InputStream is = klass.getResourceAsStream(resource)) {
			Assume.assumeNotNull(is); // ignore the test if the resource can't be found
			Files.copy(is, dest, StandardCopyOption.REPLACE_EXISTING);
		}
	}
}
