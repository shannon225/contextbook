package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;

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
		// TODO: none of the peptides in this file are found in the FASTA!
		final Path tsv = getResourceAsTempFile(tmpDir, getName(), ".tsv", "/edu/washington/gs/maccoss/encyclopedia/testdata/OpenSWATH_SM4_iRT_AssayLibrary.csv");

		final LibraryInterface library = OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(tsv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		try {
			EncyclopediaTestUtils.assertValidDlib(library); // asserts that the resulting file has DLIB extension

			// update if you change the test resource
			assertEquals("TODO: Wrong number of entries", -1, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}
}