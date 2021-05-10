package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;

public class SpectronautCSVToLibraryConverterTest extends AbstractFileConverterTest {
	private static final String NAME = "SpectronautCSVToLibraryConverterTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test
	public void testParseMods() {
		String seq=SpectronautCSVToLibraryConverter.parseMods("_YNDLEKHVCSIK_");
		assertEquals("YNDLEKHVCSIK", seq);
		seq=SpectronautCSVToLibraryConverter.parseMods("_YNDLEKHVC[Carbamidomethyl (C)]SIK_");
		assertEquals("YNDLEKHVC[+57.021464]SIK", seq);
		seq=SpectronautCSVToLibraryConverter.parseMods("_[Acetyl (Protein N-term)]YNDLEKHVCSIK_");
		assertEquals("Y[+42.010565]NDLEKHVCSIK", seq);
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNull() throws Exception {
		SpectronautCSVToLibraryConverter.convertFromSpectronautCSV(null, getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexisting() throws Exception {
		final Path csv = Files.createTempFile(tmpDir, NAME, ".csv");
		Files.delete(csv);

		SpectronautCSVToLibraryConverter.convertFromSpectronautCSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test
	public void testConvertEmptyFile() throws Exception {
		final Path csv = Files.createTempFile(tmpDir, NAME, ".csv");

		final LibraryFile library = SpectronautCSVToLibraryConverter.convertFromSpectronautCSV(csv.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		library.openFile();
		try {
			EncyclopediaTestUtils.assertValidDlib(library);

			assertEquals("Wrong number of entries", 0, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}
}
