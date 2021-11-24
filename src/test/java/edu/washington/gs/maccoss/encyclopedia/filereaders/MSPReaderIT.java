package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;

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
		// Note that 5 out of 43 entries are not found in the FASTA, even though we generated
		// it from the MSP, probably because of a defect in the generating script.
		final Path msp = getResourceAsTempFile(tmpDir, getName(), ".msp", "/edu/washington/gs/maccoss/encyclopedia/testdata/NIST_format.shuffled.abridged.msp");
		final Path fasta = EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/NIST_format.shuffled.abridged.fasta", tmpDir, NAME, ".fasta");

		MSPReader.convertMSP(msp.toFile(), fasta.toFile(), out.toFile(), SearchParameterParser.getDefaultParametersObject());

		final LibraryFile library = new LibraryFile();
		library.openFile(out.toFile());
		try {
			EncyclopediaTestUtils.assertValidDlib(library);

			// update if you change the test resource
			assertEquals("Wrong number of entries", 43 - 5, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	@Test
	public void testConvertSpTxtToLibrary() throws Exception {
		// Note that 820 out of 1319 entries are not found in the FASTA from below!
		final Path sptxt = getResourceAsTempFile(tmpDir, getName(), ".sptxt", "/edu/washington/gs/maccoss/encyclopedia/testdata/A03_human_TripleTOF.sptxt");
		final Path fasta = EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tmpDir, NAME, ".fasta");

		MSPReader.convertMSP(sptxt.toFile(), fasta.toFile(), out.toFile(), SearchParameterParser.getDefaultParametersObject());

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
}