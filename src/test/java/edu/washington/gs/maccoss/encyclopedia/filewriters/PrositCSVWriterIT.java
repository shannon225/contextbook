package edu.washington.gs.maccoss.encyclopedia.filewriters;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.Assert.assertTrue;

public final class PrositCSVWriterIT {
	public static final String NAME = "PrositCSVWriterIT";

	private Path tmpDir;
	private Path csv;

	@Before
	public void setUp() throws Exception {
		tmpDir = Files.createTempDirectory(NAME);
		FileUtils.forceDeleteOnExit(tmpDir.toFile());

		csv = Files.createTempFile(tmpDir, NAME, ".prosit.csv");
		Files.delete(csv);
	}

	@After
	public void tearDown() throws Exception {
		if (null != tmpDir) {
			// recursively delete the whole directory
			FileUtils.deleteQuietly(tmpDir.toFile());
			tmpDir = null;
			csv = null;
		}
	}

	//TODO: test fasta -> prosit

	//TODO: test elib -> prosit (incl. error states)
}