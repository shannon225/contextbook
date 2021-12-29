package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

public class TraMLToLibraryConverterTest extends AbstractFileConverterTest {
	public static final String NAME = "TraMLToLibraryConverterTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNull() throws Exception {
		TraMLSAXToLibraryConverter.convertTraML(null, getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexisting() throws Exception {
		final Path traml = Files.createTempFile(tmpDir, NAME, ".traml");
		Files.delete(traml);

		TraMLSAXToLibraryConverter.convertTraML(traml.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertEmptyFile() throws Exception {
		final Path traml = Files.createTempFile(tmpDir, NAME, ".traml");

		TraMLSAXToLibraryConverter.convertTraML(traml.toFile(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}

	@Test
	public void testConvertTraMLToLibrary() throws Exception {
		//InputStream is=getClass().getResourceAsStream("/edu/washington/gs/maccoss/encyclopedia/filereaders/ToyExample1.TraML");
		final Path traml = EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/filereaders/ToyExample1.TraML", tmpDir, NAME, ".TraML");
		final Path dest = Files.createTempFile(tmpDir, NAME, ".dlib");
		final LibraryFile library = TraMLSAXToLibraryConverter.convertTraML(traml.toFile(), null, dest.toFile(), SearchParameterParser.getDefaultParametersObject());

		// The library will be closed after conversion, so we must reopen it
		library.openFile();
		try {
			EncyclopediaTestUtils.assertValidDlib(library);

			// update if you change the test resource
			ArrayList<LibraryEntry> allEntries = library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable());
			assertEquals("Wrong number of entries", 1, allEntries.size());
			LibraryEntry entry=allEntries.get(0);
			assertEquals("wrong peptide sequence", "A[+143.058243]DTHFLLNIYDQLR", entry.getPeptideModSeq());
			assertEquals("wrong fragment number", 2, entry.getMassArray().length);
			assertEquals("wrong fragment sum", 1583.77, General.sum(entry.getMassArray()), 0.01);
			assertEquals("wrong fragment intensity sum", Float.MIN_VALUE+Float.MIN_VALUE, General.sum(entry.getIntensityArray()), 0.01f);
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}
}