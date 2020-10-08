package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Assume;
import org.junit.Before;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

import static org.junit.Assert.assertTrue;

public final class PrositCSVWriterIT {
	public static final String NAME = "PrositCSVWriterIT";

	private Path tmpDir;
	private Path csv;

	@Before
	public void setUp() throws Exception {
		tmpDir = Files.createTempDirectory(NAME);
		FileUtils.forceDeleteOnExit(tmpDir.toFile());

		csv = Files.createTempFile(tmpDir, NAME, ".prosit.csv");
		Files.delete(csv);
	}

	@After
	public void tearDown() throws Exception {
		if (null != tmpDir) {
			// recursively delete the whole directory
			FileUtils.deleteQuietly(tmpDir.toFile());
			tmpDir = null;
			csv = null;
		}
	}

	@Test
	public void testFastaToPrositCSV() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");
		try (InputStream fastaStream = getClass().getResourceAsStream("/ecoli-190209-contam_correctNL.fasta")) {
			Assume.assumeNotNull(fastaStream); // ignore the test if the resource can't be found
			Files.copy(fastaStream, fasta, StandardCopyOption.REPLACE_EXISTING);
		}

		PrositCSVWriterTest.runFastaToCsv(fasta, csv);

		assertTrue("CSV output didn't exist!", Files.exists(csv));
		assertTrue("CSV output was empty!", Files.size(csv) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(csv.toFile()))) {
			assertTrue("Got only header line in output CSV!", r.lines().count() > 1);
		}
	}

	@Test
	public void testElibToPrositCSV() throws Exception {
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");
		// TODO: use an actual resource name instead of a made-up one
		try (InputStream elibStream = getClass().getResourceAsStream("/edu/washington/gs/maccoss/encyclopedia/testdata/simple.elib")) {
			Assume.assumeNotNull(elibStream); // ignore the test if the resource can't be found
			Files.copy(elibStream, elib, StandardCopyOption.REPLACE_EXISTING);
		}

		final LibraryFile libraryFile = new LibraryFile();
		try {
			libraryFile.openFile(elib.toFile());
			PrositCSVWriterTest.runElibToCsv(libraryFile, csv);
		} finally {
			libraryFile.close();
		}

		assertTrue("CSV output didn't exist!", Files.exists(csv));
		assertTrue("CSV output was empty!", Files.size(csv) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(csv.toFile()))) {
			assertTrue("Got only header line in output CSV!", r.lines().count() > 1);
		}
	}
}