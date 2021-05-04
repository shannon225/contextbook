package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.assertTrue;

public class EndToEndIT {

	File diaFile;
	File fastaFile;
	LibraryInterface libraryInterface;
	LibraryScoringFactory libraryScoringFactory;
	SearchParameters parameters;
	Path tempDir;

	@Before
	public void setUp() throws Exception {
		parameters = SearchParameterParser.getDefaultParametersObject();
		String name = "EndToEnd";
		tempDir = Files.createTempDirectory(name);
		FileUtils.forceDeleteOnExit(tempDir.toFile());

		libraryInterface = BlibToLibraryConverter.getFile(getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/pan_human_library_600to603.dlib", tempDir, name, ".elib").toFile());
		diaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/bcs_2020jan16_600to603_hela_clib.dia", tempDir, name, ".dia").toFile();
		fastaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/pan_human_library_600to603.fasta", tempDir, name, ".fasta").toFile();
		libraryScoringFactory = new EncyclopediaOneScoringFactory(parameters);
	}

	@After
	public void tearDown() throws Exception {
		parameters = null;

		if (null != libraryInterface) {
			libraryInterface.getSourceFiles().forEach(p -> FileUtils.deleteQuietly(p.toFile()));
			libraryInterface = null;
		}
		if (null != diaFile) {
			FileUtils.deleteQuietly(diaFile);
			diaFile = null;
		}
		if (null != fastaFile) {
			FileUtils.deleteQuietly(fastaFile);
			fastaFile = null;
		}
		if (null != tempDir) {
			FileUtils.deleteDirectory(tempDir.toFile());
			tempDir = null;
		}

		libraryScoringFactory = null;
	}

	@Test
	public void testWholePipeline() throws Exception {
		Encyclopedia.runSearch(new EmptyProgressIndicator(),new EncyclopediaJobData(diaFile,fastaFile,libraryInterface,libraryScoringFactory));
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		File outputFile = FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib");

		
	}
}
