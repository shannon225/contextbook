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
		// Note that 81 out of 96 entries are not found in the FASTA from below!
		final Path csv = getResourceAsTempFile(tmpDir, getName(), ".csv", "/edu/washington/gs/maccoss/encyclopedia/testdata/human.abridged.spectronaut");

		final LibraryFile library = SpectronautCSVToLibraryConverter.convertFromSpectronautCSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());

		// The library will be closed after conversion, so we must reopen it
		library.openFile();
		try {
			EncyclopediaTestUtils.assertValidDlib(library); // asserts that the resulting file has DLIB extension

			// update if you change the test resource -- grep -Po '_([A-Z]+|\[[^\]]+\])+_' src/test/resources/edu/washington/gs/maccoss/encyclopedia/testdata/coronavirus.abridged.spectronaut | sort | uniq | wc -l
			assertEquals("Wrong number of entries", 96 - 81, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tmpDir, NAME, ".fasta");
	}
}