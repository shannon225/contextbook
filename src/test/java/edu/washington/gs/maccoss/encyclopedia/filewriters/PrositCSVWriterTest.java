package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public final class PrositCSVWriterTest extends AbstractFileConverterTest {
	public static final String NAME = "PrositCSVWriterTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return null;
	}

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

	@Test
	public void testNonexistFasta() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");
		Files.delete(fasta);

		runFastaToCsv(fasta, out);

		// This is treated the same as an empty file

		assertTrue("CSV output didn't exist!", Files.exists(out));
		assertTrue("CSV output was empty!", Files.size(out) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertEquals("Got more than the expected header line in CSV output",
					1, r.lines().count()
			);
		}
	}

	@Test
	public void testEmptyFasta() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");

		runFastaToCsv(fasta, out);

		assertTrue("CSV output didn't exist!", Files.exists(out));
		assertTrue("CSV output was empty!", Files.size(out) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertEquals("Got more than the expected header line in CSV output",
					1, r.lines().count()
			);
		}
	}

	//TODO: test elib -> prosit (incl. error states)

	static void runElibToCsv(LibraryFile libraryFile, Path csv) throws Exception {
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
		runElibToCsv(null, out);
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
			runElibToCsv(libraryFile, out);
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
		} finally {
			// clean up temp file
			libraryFile.close();
		}

		assertTrue("CSV output didn't exist!", Files.exists(csv));
		assertTrue("CSV output was empty!", Files.size(csv) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(csv.toFile()))) {
			assertEquals("Got more than the expected header line in CSV output",
					1, r.lines().count()
			);
		}
	}
}