package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;

public class MS2PIPReaderIT extends AbstractFileConverterTest {
	public static final String NAME = "MS2PIPReaderIT";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test
	public void testConvertMspToLibrary() throws Exception {
		// TODO: use an actual resource names instead of made-up ones
		final Path peps = getResourceAsTempFile(tmpDir, getName(), ".peps", "/edu/washington/gs/maccoss/encyclopedia/testdata/ms2pip/simple.peps");
		final Path csv = getResourceAsTempFile(tmpDir, getName(), ".peps", "/edu/washington/gs/maccoss/encyclopedia/testdata/ms2pip/simple.csv");

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