package edu.washington.gs.maccoss.encyclopedia.filereaders;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.Test;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class OpenSwathTSVToLibraryConverterTest extends AbstractFileConverterTest {
	private static final String NAME = "OpenSwathTSVToLibraryConverterTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test
	public void testParseUnimodMods() {
		assertEquals("P[+42.010565]EPC[+57.021464]PEPM[+15.994915]PEPR[-0.984016]", OpenSwathTSVToLibraryConverter.parseMods(".(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2)"));
	}

	@Test
	public void testParseMods() {
		// The comments in the impl list several different example cases:
		//TODO: verify if the examples here cover the observed formats in real-world files

		// TPP:    n[43]PEPC[160]PEPM[147]PEPRc[16]
		String seq=OpenSwathTSVToLibraryConverter.parseMods("n[43]PEPC[160]PEPM[147]PEPRc[16]");
		assertEquals("[43]PEPC[160]PEPM[147]PEPR[16]", seq);

		// Unimod: .(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2) (but no mods)
		seq=OpenSwathTSVToLibraryConverter.parseMods(".PEPCPEPMPEPR.");
		assertEquals("PEPCPEPMPEPR", seq);

		// TPP:    n[43]PEPC[160]PEPM[147]PEPRc[16] (but no mods)
		seq=OpenSwathTSVToLibraryConverter.parseMods("nPEPCPEPMPEPRc");
		assertEquals("PEPCPEPMPEPR", seq);

		// Fall-through case
		seq=OpenSwathTSVToLibraryConverter.parseMods("[43]PEPC[160]PEPM[147]PEPR[16]");
		assertEquals("[43]PEPC[160]PEPM[147]PEPR[16]", seq);

		seq=OpenSwathTSVToLibraryConverter.parseMods("PEPCPEPMPEPR");
		assertEquals("PEPCPEPMPEPR", seq);
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNull() throws Exception {
		OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(null, getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexisting() throws Exception {
		final Path csv = Files.createTempFile(tmpDir, NAME, ".csv");
		Files.delete(csv);

		OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test
	public void testConvertEmptyFile() throws Exception {
		final Path csv = Files.createTempFile(tmpDir, NAME, ".csv");

		final LibraryInterface library = OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
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
