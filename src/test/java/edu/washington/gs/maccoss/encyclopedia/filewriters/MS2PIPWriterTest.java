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

	@Test
	public void testGetPeptideModSeq() {
		String seq= MS2PIPWriter.getPeptideModSeq("APEPTIDEK", "-");
		assertEquals("APEPTIDEK", seq);

		seq= MS2PIPWriter.getPeptideModSeq("APEPTIDEK", "0|Deamidated");
		assertEquals("A[0.9840155826305974]PEPTIDEK", seq);

		seq= MS2PIPWriter.getPeptideModSeq("APEPTIDEK", "1|Deamidated");
		assertEquals("A[0.9840155826305974]PEPTIDEK", seq);

		seq= MS2PIPWriter.getPeptideModSeq("APEPTIDEK", "9|Oxidation");
		assertEquals("APEPTIDEK[15.994915]", seq);

		seq= MS2PIPWriter.getPeptideModSeq("APEPTIDEK", "1|Deamidated|9|Oxidation");
		assertEquals("A[0.9840155826305974]PEPTIDEK[15.994915]", seq);
	}

	@Test(expected = NullPointerException.class)
	public void testCheckPEPRECNameNullFile() throws Exception {
		MS2PIPWriter.checkPEPRECName(
				null,
				null,
				DigestionEnzyme.getEnzyme("trypsin")
		);
	}

	@Test(expected = NullPointerException.class)
	public void testCheckPEPRECNameNullFile2() throws Exception {
		MS2PIPWriter.checkPEPRECName(
				"",
				null,
				DigestionEnzyme.getEnzyme("trypsin")
		);
	}

	@Test
	public void testCheckPEPRECNameNullFile3() throws Exception {
		final String name = "test.peprec";
		assertEquals(name,
				MS2PIPWriter.checkPEPRECName(
						name,
						null,
						DigestionEnzyme.getEnzyme("trypsin")
				)
		);
	}

	@Test
	public void testCheckPEPRECName() throws Exception {
		final Path fasta = Files.createTempFile(tmpDir, NAME, ".fasta");
		Files.delete(fasta);

		assertEquals(fasta.toString() + ".trypsin.peprec",
				MS2PIPWriter.checkPEPRECName(
						null,
						fasta.toFile(),
						DigestionEnzyme.getEnzyme("trypsin")
				)
		);
	}

	static void runFastaToPeprec(Path fasta, Path peprec) throws FileNotFoundException {
		runFastaToPeprec(fasta.toFile(), peprec.toFile());
	}

	static void runFastaToPeprec(File fasta, File peprec) throws FileNotFoundException {
		MS2PIPWriter.writeMS2PIP(
				peprec.getAbsolutePath(),
				fasta,
				DigestionEnzyme.getEnzyme("trypsin"),
				(byte) 2,
				(byte) 2,
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

		//TODO: this assertion seems to not match the PEPREC format
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

		//TODO: this assertion seems to not match the PEPREC format
		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertEquals("Got more than the expected header line in PEPREC output",
					1, r.lines().count()
			);
		}
	}

	static void runElibToPeprec(LibraryFile libraryFile, Path peprec) throws Exception {
		MS2PIPWriter.writeMS2PIP(
				peprec.toString(),
				libraryFile,
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

		//TODO: this assertion seems to not match the PEPREC format
		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertEquals("Got more than the expected header line in PEPREC output",
					1, r.lines().count()
			);
		}
	}
}
