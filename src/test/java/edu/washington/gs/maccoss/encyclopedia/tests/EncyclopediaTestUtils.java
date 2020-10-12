package edu.washington.gs.maccoss.encyclopedia.tests;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.junit.Assume;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

import static org.junit.Assert.*;

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

	public static void cleanupLibrary(LibraryInterface library) {
		if (library instanceof LibraryFile) {
			((LibraryFile) library).close();
			FileUtils.deleteQuietly(((LibraryFile) library).getFile());
		}
	}

	public static void assertValidDlib(LibraryInterface library) {
		assertNotNull(library);

		final File file = ((LibraryFile) library).getFile();
		assertEquals(LibraryFile.DLIB, "." + FilenameUtils.getExtension(file.getName()));

		assertTrue(Files.exists(file.toPath()));
	}
}
