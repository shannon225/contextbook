package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;

public class LibraryToBlibConverterIT extends AbstractFileConverterTest {
	public static final String NAME = "TraMLToLibraryConverterIT";

	@Override
	protected String getName() {
		return NAME;
	}

	@Override
	protected String getOutputExtension() {
		return BlibFile.BLIB;
	}

	@Test
	public void testConvertMspToLibrary() throws Exception {
		// TODO: use an actual resource name instead of a made-up one
		final Path elib = getResourceAsTempFile(tmpDir, getName(), ".elib", "/edu/washington/gs/maccoss/encyclopedia/testdata/simple.elib");

		LibraryToBlibConverter.convert(elib.toFile(), out.toFile());

		final BlibFile blib = new BlibFile();
		blib.openFile(out.toFile());
		try {
			EncyclopediaTestUtils.assertValidBlib(blib);
		} finally {
			EncyclopediaTestUtils.cleanupBlib(blib);
		}
	}
}