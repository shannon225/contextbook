package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;

public class TraMLToLibraryConverterIT extends AbstractFileConverterTest {
	public static final String NAME = "TraMLToLibraryConverterIT";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test
	public void testConvertTraMLToLibrary() throws Exception {
		// TODO: none of the peptides in this file are found in the FASTA!
		final Path traml = getResourceAsTempFile(tmpDir, getName(), ".traml", "/edu/washington/gs/maccoss/encyclopedia/testdata/OpenSWATH_SM4_iRT_AssayLibrary.TraML");

		final LibraryFile library = TraMLToLibraryConverter.convertTraML(traml.toFile(), getFasta().toFile(), out.toFile(), SearchParameterParser.getDefaultParametersObject());

		// The library will be closed after conversion, so we must reopen it
		library.openFile();
		try {
			EncyclopediaTestUtils.assertValidDlib(library);

			// update if you change the test resource
			assertEquals("TODO: Wrong number of entries", 10, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/OpenSWATH_SM4_iRT_AssayLibrary.fasta", tmpDir, NAME, ".fasta");
	}
}