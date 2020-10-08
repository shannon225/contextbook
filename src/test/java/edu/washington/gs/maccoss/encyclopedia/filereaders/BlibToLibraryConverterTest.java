package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.tools.ant.taskdefs.optional.extension.LibFileSet;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class BlibToLibraryConverterTest {
	public static final String NAME = "BlibToLibraryConverterTest";

	private Path tmpDir;
	private Path lib;

	@Before
	public void setUp() throws Exception {
		tmpDir = Files.createTempDirectory(NAME);
		FileUtils.forceDeleteOnExit(tmpDir.toFile());

		lib = Files.createTempFile(tmpDir, NAME, "dlib");
		Files.delete(lib);
		FileUtils.forceDeleteOnExit(lib.toFile());
	}

	@After
	public void tearDown() throws Exception {
		if (null != tmpDir) {
			// recursively delete the whole directory
			FileUtils.deleteQuietly(tmpDir.toFile());
			tmpDir = null;
			lib = null;
		}
	}

	@Test(expected = NullPointerException.class)
	public void getFileFromNull() throws Exception {
		BlibToLibraryConverter.getFile(null);
	}

	@Test(expected = EncyclopediaException.class)
	public void getFileFromNonexisting() throws Exception {
		FileUtils.deleteQuietly(lib.toFile());

		BlibToLibraryConverter.getFile(lib.toFile());
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
			Assert.assertNotNull("Got null library from ELIB", library);
		} finally {
			if (library instanceof LibraryFile) {
				((LibraryFile) library).close();
				FileUtils.deleteQuietly(((LibraryFile) library).getFile());
			}
		}
	}

	@Test
	public void getFileWithExistingElib() throws Exception {
		//TODO: use a real elib instead of a 0-byte file
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");
		final Path file = tmpDir.resolve(FilenameUtils.getBaseName(elib.toFile().getName()) + ".txt");

		final LibraryInterface library = BlibToLibraryConverter.getFile(file.toFile());
		try {
			Assert.assertNotNull("Got null library from ELIB", library);
		} finally {
			if (library instanceof LibraryFile) {
				((LibraryFile) library).close();
				FileUtils.deleteQuietly(((LibraryFile) library).getFile());
			}
		}
	}

	//TODO: test dlib
	//TODO: test existing dlib

	//TODO: test conversion output auto naming (dlib extension)

	//TODO: test null/nonexist/empty conversion
}