package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;

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
	public void testConvertMs2PipToLibrary() throws Exception {
		final Path peprec = getResourceAsTempFile(tmpDir, getName(), ".peprec", "/edu/washington/gs/maccoss/encyclopedia/testdata/SGS_AQUAProteins.fasta.trypsin.peprec");
		final Path csv = getResourceAsTempFile(tmpDir, getName(), ".csv", "/edu/washington/gs/maccoss/encyclopedia/testdata/SGS_AQUAProteins.fasta.trypsin.ms2pip_predictions.csv");

		final LibraryFile library = MS2PIPReader.convertMS2PIP(peprec.toFile(), csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		try {
			EncyclopediaTestUtils.assertValidDlib(library);

			// update if you change the test resource -- one less than `wc -l SGS_AQUAProteins.fasta.trypsin.peprec`
			assertEquals("Wrong number of entries", 874, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/SGS_AQUAProteins.fasta", tmpDir, NAME, ".fasta");
	}
}