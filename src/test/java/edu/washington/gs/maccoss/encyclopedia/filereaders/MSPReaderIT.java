package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Optional;

import static org.junit.Assert.assertEquals;

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
	public void testConvertMspToLibrary() throws Exception {
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

	@Test
	public void testConvertSpTxtToLibrary() throws Exception {
		// Note that 820 out of 1319 entries are not found in the FASTA from below!
		final Path sptxt = getResourceAsTempFile(tmpDir, getName(), ".sptxt", "/edu/washington/gs/maccoss/encyclopedia/testdata/A03_human_TripleTOF.sptxt");

		MSPReader.convertMSP(sptxt.toFile(), getFasta().toFile(), out.toFile(), SearchParameterParser.getDefaultParametersObject());

		final LibraryFile library = new LibraryFile();
		library.openFile(out.toFile());
		try {
			EncyclopediaTestUtils.assertValidDlib(library);

			// update if you change the test resource
			assertEquals("Wrong number of entries", 1319 - 820, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tmpDir, NAME, ".fasta");
	}
}