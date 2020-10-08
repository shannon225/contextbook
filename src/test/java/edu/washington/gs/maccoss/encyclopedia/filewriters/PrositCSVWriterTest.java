package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.*;

public final class PrositCSVWriterTest {
	public static final String NAME = "PrositCSVWriterTest";

	private Path tmpDir;

	@Before
	public void setUp() throws Exception {
		tmpDir = Files.createTempDirectory(NAME);
		FileUtils.forceDeleteOnExit(tmpDir.toFile());
	}

	@After
	public void tearDown() throws Exception {
		if (null != tmpDir) {
			FileUtils.deleteQuietly(tmpDir.toFile());
			tmpDir = null;
		}
	}

	//TODO: test fasta -> prosit (incl. error states)

	@Test(expected = NullPointerException.class)
	public void testNullFasta() throws Exception {
		PrositCSVWriter.writeCSV(null);
	}

	@Test(expected = IOException.class)
	public void testNonexistFasta() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");
		Files.delete(fasta);

		PrositCSVWriter.writeCSV(fasta.toFile());
	}

	//TODO: test elib -> prosit (incl. error states)

	@Test(expected = NullPointerException.class)
	public void testNullElib() throws Exception {
		PrositCSVWriter.writeCSV((LibraryFile) null, 50, (byte) 2, true);
	}

	@Test(expected = NullPointerException.class)
	public void testNullElibFile() throws Exception {
		PrositCSVWriter.writeCSV(new LibraryFile(), 50, (byte) 2, true);
	}

	@Test(expected = IOException.class)
	@Ignore //TODO: figure this out
	public void testNonexistElib() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");
		Files.delete(fasta);

		PrositCSVWriter.writeCSV(fasta.toFile());
	}
}