package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Optional;

public class MSPReaderIT extends AbstractFileConverterTest {
	public static final String NAME = "MSPReaderIT";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test
	public void testConvertBlibToLibrary() throws Exception {
		// TODO: use an actual resource name instead of a made-up one
		final Path msp = getResourceAsTempFile(tmpDir, getName(), ".msp", "/edu/washington/gs/maccoss/encyclopedia/testdata/simple.msp");

		MSPReader.convertMSP(msp.toFile(), getFasta().toFile(), out.toFile(), SearchParameterParser.getDefaultParametersObject());

		final LibraryFile library = new LibraryFile();
		library.openFile(out.toFile());
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