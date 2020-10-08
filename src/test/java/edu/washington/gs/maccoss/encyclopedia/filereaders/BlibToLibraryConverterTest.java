package edu.washington.gs.maccoss.encyclopedia.filereaders;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;

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

	@Test
	public void getFile() {
		//TODO
	}

	//TODO: test conversion output auto naming (dlib extension)

	//TODO: test null/nonexist/empty conversion
}