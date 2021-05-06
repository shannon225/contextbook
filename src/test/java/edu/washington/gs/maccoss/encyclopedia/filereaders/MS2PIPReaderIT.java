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
	public void testConvertMs2PipToLibrary() throws Exception {
		final Path peprec = getResourceAsTempFile(tmpDir, getName(), ".peprec", "/edu/washington/gs/maccoss/encyclopedia/testdata/SGS_AQUAProteins.fasta.trypsin.peprec");
		final Path csv = getResourceAsTempFile(tmpDir, getName(), ".csv", "/edu/washington/gs/maccoss/encyclopedia/testdata/SGS_AQUAProteins.fasta.trypsin.ms2pip_predictions.csv");

		final LibraryFile library = MS2PIPReader.convertMS2PIP(peprec.toFile(), csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
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