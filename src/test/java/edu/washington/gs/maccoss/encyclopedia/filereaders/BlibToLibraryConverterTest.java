package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
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
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");
		{ // set up elib schema
			final LibraryFile libraryFile = new LibraryFile();
			try {
				libraryFile.openFile(elib.toFile());
				libraryFile.saveFile();
			} finally {
				libraryFile.close();
			}
		}

		final LibraryInterface library = BlibToLibraryConverter.getFile(elib.toFile());
		try {
			assertNotNull("Got null library from ELIB", library);
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	@Test
	public void getFileWithExistingElib() throws Exception {
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");
		{ // set up elib schema
			final LibraryFile libraryFile = new LibraryFile();
			try {
				libraryFile.openFile(elib.toFile());
				libraryFile.saveFile();
			} finally {
				libraryFile.close();
			}
		}

		final Path file = tmpDir.resolve(FilenameUtils.getBaseName(elib.toFile().getName()) + ".txt");
		Files.write(file, new byte[0], StandardOpenOption.CREATE_NEW); // create empty file

		final LibraryInterface library = BlibToLibraryConverter.getFile(file.toFile());
		try {
			assertNotNull("Got null library from ELIB", library);
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	@Test
	public void getFileFromDlib() throws Exception {
		{ // set up elib schema
			final LibraryFile libraryFile = new LibraryFile();
			try {
				Files.write(out, new byte[0], StandardOpenOption.CREATE_NEW); // create empty file
				libraryFile.openFile(out.toFile());
				libraryFile.saveFile();
			} finally {
				libraryFile.close();
			}
		}

		final LibraryInterface library = BlibToLibraryConverter.getFile(out.toFile());
		try {
			assertNotNull("Got null library from DLIB", library);
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	@Test
	public void getFileWithExistingDlib() throws Exception {
		{ // set up elib schema
			final LibraryFile libraryFile = new LibraryFile();
			try {
				Files.write(out, new byte[0], StandardOpenOption.CREATE_NEW); // create empty file
				libraryFile.openFile(out.toFile());
				libraryFile.saveFile();
			} finally {
				libraryFile.close();
			}
		}

		final Path file = tmpDir.resolve(FilenameUtils.getBaseName(out.toFile().getName()) + ".txt");
		Files.write(file, new byte[0], StandardOpenOption.CREATE_NEW); // create empty file

		final LibraryInterface library = BlibToLibraryConverter.getFile(file.toFile());
		try {
			assertNotNull("Got null library from DLIB", library);
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	@Test
	public void testConvertEmptyBlib() throws Exception {
		final Path blib = EncyclopediaTestUtils.getResourceAsTempFile(getClass(), "/empty.blib", tmpDir, NAME, ".blib");

		final LibraryFile library = BlibToLibraryConverter.convert(blib.toFile(), Optional.empty(), getFasta().toFile(), true, SearchParameterParser.getDefaultParametersObject());
		library.openFile();
		try {
			EncyclopediaTestUtils.assertValidDlib(library); // asserts that the resulting file has DLIB extension

			assertEquals("Wrong number of entries", 0, library.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		} finally {
			EncyclopediaTestUtils.cleanupLibrary(library);
		}
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNull() throws Exception {
		BlibToLibraryConverter.convert(null, Optional.empty(), getFasta().toFile(), true, SearchParameterParser.getDefaultParametersObject());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexisting() throws Exception {
		final Path blib = Files.createTempFile(tmpDir, NAME, ".blib");
		Files.delete(blib);

		BlibToLibraryConverter.convert(blib.toFile(), Optional.empty(), getFasta().toFile(), true, SearchParameterParser.getDefaultParametersObject());
	}

	@Test
	public void testConvertEmptyFile() throws Exception {
		final Path blib = Files.createTempFile(tmpDir, NAME, ".blib");

		final LibraryFile library = BlibToLibraryConverter.convert(blib.toFile(), Optional.empty(), getFasta().toFile(), true, SearchParameterParser.getDefaultParametersObject());
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