package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertTrue;

public final class PrositCSVWriterIT extends AbstractFileConverterTest {
	@Override
	protected String getName() {
		return "PrositCSVWriterIT";
	}

	@Override
	protected String getOutputExtension() {
		return ".prosit.csv";
	}

	@Test
	public void testFastaToPrositCSV() throws Exception {
		final Path fasta = getResourceAsTempFile(this.tmpDir, getName(), ".fasta", "/ecoli-190209-contam_correctNL.fasta");

		PrositCSVWriterTest.runFastaToCsv(fasta, out);

		assertTrue("CSV output didn't exist!", Files.exists(out));
		assertTrue("CSV output was empty!", Files.size(out) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertTrue("Got only header line in output CSV!", r.lines().count() > 1);
		}
	}

	@Test
	public void testElibToPrositCSV() throws Exception {
		final Path elib = getResourceAsTempFile(tmpDir, getName(), ".elib", "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia.elib");

		final LibraryFile libraryFile = new LibraryFile();
		try {
			libraryFile.openFile(elib.toFile());
			PrositCSVWriterTest.runElibToCsv(libraryFile, out);
		} finally {
			libraryFile.close();
		}

		assertTrue("CSV output didn't exist!", Files.exists(out));
		assertTrue("CSV output was empty!", Files.size(out) > 0);

		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertTrue("Got only header line in output CSV!", r.lines().count() > 1);
		}
	}
}