package edu.washington.gs.maccoss.encyclopedia.tests;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public abstract class AbstractFileConverterTest {
	protected Path tmpDir;
	protected Path out;

	@Before
	public void setUp() throws Exception {
		tmpDir = Files.createTempDirectory(getName());
		FileUtils.forceDeleteOnExit(tmpDir.toFile());

		out = Files.createTempFile(tmpDir, getName(), getOutputExtension());
		Files.delete(out);
		FileUtils.forceDeleteOnExit(out.toFile());
	}

	@After
	public void tearDown() throws Exception {
		if (null != tmpDir) {
			// recursively delete the whole directory
			FileUtils.deleteQuietly(tmpDir.toFile());
			tmpDir = null;
			out = null;
		}
	}

	protected abstract String getName();

	protected abstract String getOutputExtension();

	protected final Path getResourceAsTempFile(Path tmpDir, String name, String suffix, String resource) throws IOException {
		return EncyclopediaTestUtils.getResourceAsTempFile(getClass(), resource, tmpDir, name, suffix);
	}
}
