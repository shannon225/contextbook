package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertEquals;

public class LibraryToBlibConverterTest extends AbstractFileConverterTest {
	public static final String NAME = "LibraryToBlibConverterTest";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return BlibFile.BLIB;
	}

	@Test(expected = NullPointerException.class)
	public void testConvertNull() throws Exception {
		LibraryToBlibConverter.convert(null, out.toFile());
	}

	@Test(expected = EncyclopediaException.class)
	public void testConvertNonexisting() throws Exception {
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");
		Files.delete(elib);

		LibraryToBlibConverter.convert(elib.toFile(), out.toFile());
	}

	@Test
	public void testConvertEmptyFile() throws Exception {
		final Path elib = Files.createTempFile(tmpDir, NAME, ".elib");

		LibraryToBlibConverter.convert(elib.toFile(), out.toFile());

		final BlibFile blib = new BlibFile();
		blib.openFile(out.toFile());
		try {
			EncyclopediaTestUtils.assertValidBlib(blib);

			EncyclopediaTestUtils.assertEmptyBlib(blib);
		} finally {
			EncyclopediaTestUtils.cleanupBlib(blib);
		}
	}
}