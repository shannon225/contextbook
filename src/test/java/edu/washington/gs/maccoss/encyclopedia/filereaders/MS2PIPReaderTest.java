package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.apache.commons.io.FileUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class MS2PIPReaderTest extends AbstractFileConverterTest {
	public static final String NAME = "MS2PIPReaderTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNullPeps() throws Exception {
		final Path csv = Files.createTempFile(tmpDir, NAME, ".csv");
		FileUtils.forceDeleteOnExit(csv.toFile());

		MS2PIPReader.convertMS2PIP(null, csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNullCsv() throws Exception {
		final Path peps = Files.createTempFile(tmpDir, NAME, ".peps");
		FileUtils.forceDeleteOnExit(peps.toFile());

		MS2PIPReader.convertMS2PIP(peps.toFile(), null, getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexistingPeps() throws Exception {
		final Path peps = Files.createTempFile(tmpDir, NAME, ".peps");
		Files.delete(peps);
		final Path csv = Files.createTempFile(tmpDir, NAME, ".csv");
		FileUtils.forceDeleteOnExit(csv.toFile());

		MS2PIPReader.convertMS2PIP(peps.toFile(), csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexistingCsv() throws Exception {
		final Path peps = Files.createTempFile(tmpDir, NAME, ".peps");
		FileUtils.forceDeleteOnExit(peps.toFile());
		final Path csv = Files.createTempFile(tmpDir, NAME, ".csv");
		Files.delete(csv);

		MS2PIPReader.convertMS2PIP(peps.toFile(), csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	// Succeeds and produces an empty library
	@Test
	public void testConvertEmptyFiles() throws Exception {
		final Path peps = Files.createTempFile(tmpDir, NAME, ".peps");
		FileUtils.forceDeleteOnExit(peps.toFile());
		final Path csv = Files.createTempFile(tmpDir, NAME, ".csv");
		FileUtils.forceDeleteOnExit(csv.toFile());

		// Succeeds
		final LibraryFile library = MS2PIPReader.convertMS2PIP(peps.toFile(), csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		try {
			EncyclopediaTestUtils.assertValidDlib(library);
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}

	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}
}