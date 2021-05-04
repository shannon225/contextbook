package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Optional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class BlibToLibraryConverterIT extends AbstractFileConverterTest {
	public static final String NAME = "BlibToLibraryConverterIT";

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
		final Path blib = getResourceAsTempFile(tmpDir, getName(), ".blib", "/edu/washington/gs/maccoss/encyclopedia/testdata/YeastDDA_subset.blib");

		final LibraryFile library = BlibToLibraryConverter.convert(blib.toFile(), Optional.empty(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());

		// The library will be closed after conversion, so we must reopen it
		library.openFile();
		try {
			EncyclopediaTestUtils.assertValidDlib(library); // asserts that the resulting file has DLIB extension

			// update if you change the test resource
			assertEquals("Wrong number of entries", 23, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		//TODO: this might be a little large
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/yeast-contam.fasta", tmpDir, NAME, ".fasta");
	}
}