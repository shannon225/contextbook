package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class TraMLToLibraryConverterTest extends AbstractFileConverterTest {
	public static final String NAME = "TraMLToLibraryConverterTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNull() throws Exception {
		TraMLToLibraryConverter.convertTraML(null, getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexisting() throws Exception {
		final Path traml = Files.createTempFile(tmpDir, NAME, ".traml");
		Files.delete(traml);

		TraMLToLibraryConverter.convertTraML(traml.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertEmptyFile() throws Exception {
		final Path traml = Files.createTempFile(tmpDir, NAME, ".traml");

		TraMLToLibraryConverter.convertTraML(traml.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}
}