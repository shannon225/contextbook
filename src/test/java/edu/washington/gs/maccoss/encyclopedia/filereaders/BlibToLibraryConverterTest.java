package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.junit.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

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

/* TODO
	@Test
	public void testConvertOutputFilename() throws Exception {
		final LibraryInterface library = BlibToLibraryConverter.convert(blib.toFile(), Optional.empty(), fasta.toFile(), SearchParameterParser.getDefaultParametersObject());
		try {
			final File file = null; //TODO
			assertEquals(LibraryFile.DLIB, "." + FilenameUtils.getExtension(file.getName()));
		} finally {
			cleanupLibrary(library);
		}
	}
*/

	//TODO: test null/nonexist/empty conversion

	static void cleanupLibrary(LibraryInterface library) {
		if (library instanceof LibraryFile) {
			((LibraryFile) library).close();
			FileUtils.deleteQuietly(((LibraryFile) library).getFile());
		}
	}
}