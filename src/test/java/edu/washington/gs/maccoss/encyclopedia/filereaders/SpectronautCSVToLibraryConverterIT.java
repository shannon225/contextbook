package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Optional;

import static org.junit.Assert.assertEquals;

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
		final Path csv = getResourceAsTempFile(tmpDir, getName(), ".csv", "/edu/washington/gs/maccoss/encyclopedia/testdata/coronavirus.abridged.spectronaut");

		final LibraryInterface library = SpectronautCSVToLibraryConverter.convertFromSpectronautCSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		try {
			EncyclopediaTestUtils.assertValidDlib(library); // asserts that the resulting file has DLIB extension

			// update if you change the test resource
			assertEquals("Wrong number of entries", 39, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		//TODO: this is surely too large
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/2020-02-27_uniprot_organism_coronavirus.fasta", tmpDir, NAME, ".fasta");
	}
}