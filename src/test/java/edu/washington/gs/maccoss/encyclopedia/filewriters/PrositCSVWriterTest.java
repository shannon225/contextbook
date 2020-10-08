package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertTrue;

public final class PrositCSVWriterTest {
	public static final String NAME = "PrositCSVWriterTest";

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

	//TODO: test fasta -> prosit (incl. error states)

	static void runFastaToCsv(Path fasta, Path csv) throws FileNotFoundException {
		runFastaToCsv(fasta.toFile(), csv.toFile());
	}

	static void runFastaToCsv(File fasta, File csv) throws FileNotFoundException {
		PrositCSVWriter.writeCSV(
				csv.getAbsolutePath(),
				fasta,
				DigestionEnzyme.getEnzyme("trypsin"),
				50,
				(byte) 2,
				(byte) 2, (byte) 2,
				1,
				new Range(400f, 1000f),
				true
		);
	}

	@Test(expected = NullPointerException.class)
	public void testNullFasta() throws Exception {
		runFastaToCsv(null, Files.createTempFile(tmpDir, NAME, ".csv").toFile());
	}

	@Test(expected = FileNotFoundException.class)
	public void testNonexistFasta() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");
		Files.delete(fasta);

		runFastaToCsv(fasta, csv);
	}

	@Test
	public void testEmptyFasta() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");

		runFastaToCsv(fasta, csv);

		assertTrue("CSV output didn't exist!", Files.exists(csv));
		assertTrue("CSV output was empty!", Files.size(csv) > 0);
	}

	//TODO: test elib -> prosit (incl. error states)

	void runElibToCsv(LibraryFile libraryFile, Path csv) throws Exception {
		PrositCSVWriter.writeCSV(
				csv.toString(),
				libraryFile,
				50,
				(byte) 2,
				true
		);
	}

	@Test(expected = NullPointerException.class)
	public void testNullElib() throws Exception {
		runElibToCsv(null, csv);
	}

	@Test(expected = NullPointerException.class)
	public void testCheckCSVNameNullElibFile() throws Exception {
		PrositCSVWriter.checkCSVName(
				null,
				null,
				DigestionEnzyme.getEnzyme("trypsin"),
				50,
				(byte) 2
		);
	}

	@Test(expected = EncyclopediaException.class)
	public void testNullElibFile() throws Exception {
		final LibraryFile libraryFile = new LibraryFile();
		try {
			runElibToCsv(libraryFile, csv);
		} finally {
			libraryFile.close();
		}
	}

	@Test
	public void testEmptyElib() throws Exception {
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");
		final Path csv = Files.createTempFile(tmpDir, NAME, ".prosit.csv");

		final LibraryFile libraryFile = new LibraryFile();
		try {
			libraryFile.openFile(elib.toFile());

			runElibToCsv(libraryFile, csv);

			assertTrue("CSV output didn't exist!", Files.exists(csv));
			assertTrue("CSV output was empty!", Files.size(csv) > 0);
		} finally {
			// clean up temp file
			libraryFile.close();
		}
	}
}