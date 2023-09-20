package edu.washington.gs.maccoss.encyclopedia.filereaders;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
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

		// Unimod w/o pre/post AAs: .(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2)
		seq=OpenSwathTSVToLibraryConverter.parseMods(".(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2)");
		assertEquals("P[+42.010565]EPC[+57.021464]PEPM[+15.994915]PEPR[-0.984016]", seq);

		// Unimod: .(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2) (but no mods)
		seq=OpenSwathTSVToLibraryConverter.parseMods(".PEPCPEPMPEPR.");
		assertEquals("PEPCPEPMPEPR", seq);

		// Unimod with pre/post AAs
		seq=OpenSwathTSVToLibraryConverter.parseMods("A.VNC(UniMod:4)IQWIC(UniMod:4)K.E");
		assertEquals("VNC[+57.021464]IQWIC[+57.021464]K", seq);

		// Unimod with pre/post AAs elided
		seq=OpenSwathTSVToLibraryConverter.parseMods("-.VNC(UniMod:4)IQWIC(UniMod:4)K.-");
		assertEquals("VNC[+57.021464]IQWIC[+57.021464]K", seq);

		// Unimod with pre/post AAs (but no mods)
		seq=OpenSwathTSVToLibraryConverter.parseMods("A.VNCIQWICK.E");
		assertEquals("VNCIQWICK", seq);

		// Unimod with pre/post AAs elided (but no mods)
		seq=OpenSwathTSVToLibraryConverter.parseMods("-.VNCIQWICK.-");
		assertEquals("VNCIQWICK", seq);

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

		final LibraryFile library = OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		library.openFile();
		try {
			EncyclopediaTestUtils.assertValidDlib(library);

			assertEquals("Wrong number of entries", 0, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	@Test
	public void testConvertWithoutGroupId() throws Exception {
		final Path csv = getResourceAsTempFile(tmpDir, "test_", ".tsv", "/edu/washington/gs/maccoss/encyclopedia/filereaders/no-group-id.tsv");

		final LibraryFile library = OpenSwathTSVToLibraryConverter.convertFromOpenSwathTSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		library.openFile();
		try {
			EncyclopediaTestUtils.assertValidDlib(library);

			// Example file contains exactly one peptide
			final ArrayList<LibraryEntry> allEntries = library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable());
			assertEquals("Wrong number of entries", 1, allEntries.size());
			assertEquals("Wrong number of transitions", 5, allEntries.iterator().next().getMassArray().length);
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}
}
