package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Optional;

public class SpectronautCSVToLibraryConverterIT extends AbstractFileConverterTest {
	public static final String NAME = "SpectronautCSVToLibraryConverterIT";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test
	public void testConvertSpectronautCSVToLibrary() throws Exception {
		// TODO: use an actual resource name instead of a made-up one
		final Path csv = getResourceAsTempFile(tmpDir, getName(), ".csv", "/edu/washington/gs/maccoss/encyclopedia/testdata/simple.csv");

		final LibraryInterface library = SpectronautCSVToLibraryConverter.convertFromSpectronautCSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		try {
			EncyclopediaTestUtils.assertValidDlib(library); // asserts that the resulting file has DLIB extension
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}
}