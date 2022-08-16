package edu.washington.gs.maccoss.encyclopedia.tests;

import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibFile;
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
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import static org.junit.Assert.*;

public abstract class EncyclopediaTestUtils {
	private EncyclopediaTestUtils() throws IllegalAccessError {
		throw new IllegalAccessError();
	}

	public static Path getEmptyTempFile(Path tmpDir, String name, String suffix) throws IOException {
		final Path dest;
		if (tmpDir==null) {
			dest=Files.createTempFile(name, suffix);
		} else {
			dest=Files.createTempFile(tmpDir, name, suffix);
		}
		FileUtils.forceDeleteOnExit(dest.toFile());
		return dest;
	}

	public static Path getResourceAsTempFile(Class<?> klass, String resource, Path tmpDir, String name, String suffix) throws IOException {
		final Path dest;
		if (tmpDir==null) {
			dest=Files.createTempFile(name, suffix);
		} else {
			dest=Files.createTempFile(tmpDir, name, suffix);
		}
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

	public static void cleanupBlib(BlibFile blib) {
		blib.close();
		FileUtils.deleteQuietly(blib.getUserFile());
	}

	public static void assertValidDlib(LibraryInterface library) {
		assertNotNull(library);

		final File file = ((LibraryFile) library).getFile();
		assertEquals(LibraryFile.DLIB, "." + FilenameUtils.getExtension(file.getName()));

		assertTrue(Files.exists(file.toPath()));
	}

	public static void assertValidBlib(BlibFile blib) {
		assertNotNull(blib);

		final File file = blib.getUserFile();
		assertEquals(BlibFile.BLIB, "." + FilenameUtils.getExtension(file.getName()));

		assertTrue(Files.exists(file.toPath()));
	}

	public static void assertEmptyBlib(BlibFile blib) throws IOException, SQLException {
		assertBlibEntries(blib, 0);
	}

	public static void assertBlibEntries(BlibFile blib, int expected) throws IOException, SQLException {
		try (Connection c = blib.getConnection(blib.getUserFile())) {
			try (PreparedStatement ps = c.prepareStatement("SELECT count() FROM RefSpectra;")) {
				try (ResultSet rs = ps.executeQuery()) {
					assertTrue(rs.next());
					final int count = rs.getInt(1);
					assertEquals("Unexpected count of entries in BLIB! Expected " + expected + " got " + count, expected, count);
				}
			}
		}
	}
}
