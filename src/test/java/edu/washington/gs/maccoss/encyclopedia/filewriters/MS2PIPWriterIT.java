package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

public class MS2PIPWriterIT extends AbstractFileConverterTest {
	private static final String NAME = "MS2PIPWriterIT";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return ".peprec";
	}

	@Test
	public void testSimpleFasta() throws Exception {
		final Path fasta = getResourceAsTempFile(this.tmpDir, getName(), ".fasta", "/ecoli-190209-contam_correctNL.fasta");

		MS2PIPWriterTest.runFastaToPeprec(fasta, out);

		assertTrue("PEPREC output didn't exist!", Files.exists(out));
		assertTrue("PEPREC output was empty!", Files.size(out) > 0);

		assertValidPeprec(out);
	}

	@Test
	public void testSimpleElib() throws Exception {
		// TODO: use an actual resource name instead of a made-up one
		final Path elib = getResourceAsTempFile(tmpDir, getName(), ".elib", "/edu/washington/gs/maccoss/encyclopedia/testdata/simple.elib");

		final LibraryFile libraryFile = new LibraryFile();
		try {
			libraryFile.openFile(elib.toFile());

			MS2PIPWriterTest.runElibToPeprec(libraryFile, out);
		} finally {
			// clean up temp file
			libraryFile.close();
		}

		assertTrue("PEPREC output didn't exist!", Files.exists(out));
		assertTrue("PEPREC output was empty!", Files.size(out) > 0);

		assertValidPeprec(out);
	}

	static void assertValidPeprec(Path out) throws IOException {
		try (BufferedReader r = new BufferedReader(new FileReader(out.toFile()))) {
			assertNotNull("Didn't get any lines from PEPREC output", r.readLine());
			assertNotNull("Got only a header line in PEPREC output", r.readLine());
		}
	}
}
