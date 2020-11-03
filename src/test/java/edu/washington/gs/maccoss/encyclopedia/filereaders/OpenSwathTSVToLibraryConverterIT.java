package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;

public class OpenSwathTSVToLibraryConverterIT extends AbstractFileConverterTest {
	public static final String NAME = "OpenSwathTSVToLibraryConverterIT";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test
	public void testConvertOpenSwathTSVToLibrary() throws Exception {
		// TODO: use an actual resource name instead of a made-up one
		final Path tsv = getResourceAsTempFile(tmpDir, getName(), ".tsv", "/edu/washington/gs/maccoss/encyclopedia/testdata/simple.openswath.tsv");

		final LibraryInterface library = OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(tsv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
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