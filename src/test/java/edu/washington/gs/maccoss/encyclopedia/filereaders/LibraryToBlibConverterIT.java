package edu.washington.gs.maccoss.encyclopedia.filereaders;

import edu.washington.gs.maccoss.encyclopedia.tests.AbstractFileConverterTest;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import org.junit.Test;

import java.nio.file.Path;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;

import static org.junit.Assert.assertEquals;

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
	public void testConvertElibToLibrary() throws Exception {
		final Path elib = getResourceAsTempFile(tmpDir, getName(), ".elib", "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia.elib");

		LibraryToBlibConverter.convert(elib.toFile(), out.toFile());

		assertValidWithCount(407); // update if you change the test resource -- SELECT count() FROM Entries;
	}

	@Test
	public void testConvertDlibToLibrary() throws Exception {
		final Path elib = getResourceAsTempFile(tmpDir, getName(), ".dlib", "/edu/washington/gs/maccoss/encyclopedia/testdata/truncated_pan_human_library.dlib");

		LibraryToBlibConverter.convert(elib.toFile(), out.toFile());

		assertValidWithCount(4669); // update if you change the test resource -- SELECT count() FROM Entries;
	}

	private void assertValidWithCount(int count) throws Exception {
		final BlibFile blib = new BlibFile();
		blib.openFile(out.toFile());
		try {
			EncyclopediaTestUtils.assertValidBlib(blib);

			try (Connection c = blib.getConnection(out.toFile())) {
				try (PreparedStatement s = c.prepareStatement("SELECT count() from RefSpectra;")) {
					try (ResultSet rs = s.executeQuery()) {
						assertEquals("Wrong number of entries", count, rs.getInt(1));
					}
				}
			}
		} finally {
			EncyclopediaTestUtils.cleanupBlib(blib);
		}
	}
}