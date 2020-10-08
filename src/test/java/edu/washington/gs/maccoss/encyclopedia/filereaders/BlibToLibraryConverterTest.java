package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;

import static org.junit.Assert.*;

public class BlibToLibraryConverterTest extends AbstractFileConverterTest {
	public static final String NAME = "BlibToLibraryConverterTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return LibraryFile.DLIB;
	}

	@Test(expected = NullPointerException.class)
	public void getFileFromNull() throws Exception {
		BlibToLibraryConverter.getFile(null);
	}

	@Test(expected = EncyclopediaException.class)
	public void getFileFromNonexisting() throws Exception {
		FileUtils.deleteQuietly(out.toFile());

		BlibToLibraryConverter.getFile(out.toFile());
	}

	@Test(expected = EncyclopediaException.class)
	public void getFileWithoutExisting() throws Exception {
		final Path file = Files.createTempFile(tmpDir, NAME, ".txt");

		BlibToLibraryConverter.getFile(file.toFile());
	}

	@Test
	public void getFileFromElib() throws Exception {
		//TODO: use a real elib instead of a 0-byte file
		final LibraryInterface library = BlibToLibraryConverter.getFile(Files.createTempFile(tmpDir, NAME, ".elib").toFile());
		try {
			assertNotNull("Got null library from ELIB", library);
		} finally {
			cleanupLibrary(library);
		}
	}

	@Test
	public void getFileWithExistingElib() throws Exception {
		//TODO: use a real elib instead of a 0-byte file
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");
		final Path file = tmpDir.resolve(FilenameUtils.getBaseName(elib.toFile().getName()) + ".txt");

		final LibraryInterface library = BlibToLibraryConverter.getFile(file.toFile());
		try {
			assertNotNull("Got null library from ELIB", library);
		} finally {
			cleanupLibrary(library);
		}
	}

	@Test
	public void getFileFromDlib() throws Exception {
		//TODO: use a real dlib instead of a 0-byte file
		final LibraryInterface library = BlibToLibraryConverter.getFile(out.toFile());
		try {
			assertNotNull("Got null library from DLIB", library);
		} finally {
			cleanupLibrary(library);
		}
	}

	@Test
	public void getFileWithExistingDlib() throws Exception {
		//TODO: use a real dlib instead of a 0-byte file
		final Path file = tmpDir.resolve(FilenameUtils.getBaseName(out.toFile().getName()) + ".txt");

		final LibraryInterface library = BlibToLibraryConverter.getFile(file.toFile());
		try {
			assertNotNull("Got null library from DLIB", library);
		} finally {
			cleanupLibrary(library);
		}
	}

	@Test
	public void testConvertEmptyBlib() throws Exception {
		final Path blib = EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/empty.blib", tmpDir, NAME, ".blib");

		final LibraryInterface library = BlibToLibraryConverter.convert(blib.toFile(), Optional.empty(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		try {
			assertValidDlib(library); // asserts that the resulting file has DLIB extension
		} finally {
			cleanupLibrary(library);
		}
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNull() throws Exception {
		BlibToLibraryConverter.convert(null, Optional.empty(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexisting() throws Exception {
		final Path blib = Files.createTempFile(tmpDir, NAME, ".blib");
		Files.delete(blib);

		BlibToLibraryConverter.convert(blib.toFile(), Optional.empty(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
	}

	@Test
	public void testConvertEmptyFile() throws Exception {
		final Path blib = Files.createTempFile(tmpDir, NAME, ".blib");

		final LibraryInterface library = BlibToLibraryConverter.convert(blib.toFile(), Optional.empty(), getFasta().toFile(), SearchParameterParser.getDefaultParametersObject());
		try {
			assertValidDlib(library);
		} finally {
			cleanupLibrary(library);
		}
	}

	Path getFasta() throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/ecoli-190209-contam_correctNL.fasta", tmpDir, NAME, ".fasta");
	}

	static void assertValidDlib(LibraryInterface library) {
		assertNotNull(library);

		final File file = ((LibraryFile) library).getFile();
		assertEquals(LibraryFile.DLIB, "." + FilenameUtils.getExtension(file.getName()));

		assertTrue(Files.exists(file.toPath()));
	}

	static void cleanupLibrary(LibraryInterface library) {
		if (library instanceof LibraryFile) {
			((LibraryFile) library).close();
			FileUtils.deleteQuietly(((LibraryFile) library).getFile());
		}
	}
}