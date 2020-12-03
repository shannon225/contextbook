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

public class MS2PIPWriterTest extends AbstractFileConverterTest {
	private static final String NAME = "MS2PIPWriterTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return ".peprec";
	}

	//TODO: test peptidemodseq parsing

	@Test(expected = NullPointerException.class)
	public void testCheckPEPRECNameNullElibFile() throws Exception {
		MS2PIPWriter.checkPEPRECName(
				null,
				null,
				DigestionEnzyme.getEnzyme("trypsin"),
				50,
				(byte) 2
		);
	}

	static void runFastaToPeprec(Path fasta, Path peprec) throws FileNotFoundException {
		runFastaToPeprec(fasta.toFile(), peprec.toFile());
	}

	static void runFastaToPeprec(File fasta, File peprec) throws FileNotFoundException {
		MS2PIPWriter.writeCSV(
				peprec.getAbsolutePath(),
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
		runFastaToPeprec(null, out.toFile());
	}

	@Test
	public void testNonexistFasta() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");
		Files.delete(fasta);

		runFastaToPeprec(fasta, out);

		// This is treated the same as an empty file

		assertTrue("PEPREC output didn't exist!", Files.exists(out));
		assertTrue("PEPREC output was empty!", Files.size(out) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertEquals("Got more than the expected header line in PEPREC output",
					1, r.lines().count()
			);
		}
	}

	@Test
	public void testEmptyFasta() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");

		runFastaToPeprec(fasta, out);

		assertTrue("PEPREC output didn't exist!", Files.exists(out));
		assertTrue("PEPREC output was empty!", Files.size(out) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertEquals("Got more than the expected header line in PEPREC output",
					1, r.lines().count()
			);
		}
	}

	static void runElibToPeprec(LibraryFile libraryFile, Path peprec) throws Exception {
		MS2PIPWriter.writeCSV(
				peprec.toString(),
				libraryFile,
				50,
				(byte) 2,
				true
		);
	}

	@Test(expected = NullPointerException.class)
	public void testNullElib() throws Exception {
		runElibToPeprec(null, out);
	}

	@Test(expected = EncyclopediaException.class)
	public void testNullElibFile() throws Exception {
		final LibraryFile libraryFile = new LibraryFile();
		try {
			runElibToPeprec(libraryFile, out);
		} finally {
			libraryFile.close();
		}
	}

	@Test
	public void testEmptyElib() throws Exception {
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");

		final LibraryFile libraryFile = new LibraryFile();
		try {
			libraryFile.openFile(elib.toFile());

			runElibToPeprec(libraryFile, out);
		} finally {
			// clean up temp file
			libraryFile.close();
		}

		assertTrue("PEPREC output didn't exist!", Files.exists(out));
		assertTrue("PEPREC output was empty!", Files.size(out) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertEquals("Got more than the expected header line in PEPREC output",
					1, r.lines().count()
			);
		}
	}
}
