package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorVersion;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ScoringBreadthType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchToELIBJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import gnu.trove.map.hash.TCharDoubleHashMap;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.LoggerFactory;

import java.awt.*;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.concurrent.ForkJoinPool;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.*;

public class EndToEndIT {

	File diaFile;
	File diaFile2;
	File diaFile3;

	File fastaFile;
	File libraryFile;
	LibraryFile libraryInterface;
	LibraryScoringFactory libraryScoringFactory;
	SearchParameters parameters;
	Path tempDir;

	static Range STANDARD_RANGE = new Range(592.5840338877389,604.3740813086648);
	static int MAX_POSSIBLE_PEPTIDES = 4669;
	static int MAX_POSSIBLE_PROTEIN_GROUPS = 6676;

	@Before
	public void setUp() throws Exception {
		if (GraphicsEnvironment.isHeadless() && "1.8".equals(System.getProperty("java.specification.version"))) {
			LoggerFactory.getLogger(EndToEndIT.class)
					.info("Disabling assistive technologies to avoid errors in headless build on Java 8!");

			// On JDK 8 running headless we can encounter problems if this is set by the system
			// installation; by overwriting it we avoid errors with JFreeChart. This is not an
			// issue on later releases that will better respect the headless flag.
			// See https://stackoverflow.com/a/59397731/115714
			System.setProperty("javax.accessibility.assistive_technologies", "java.lang.Object");
		}

		SearchParameters parameters = SearchParameterParser.getDefaultParametersObject();
		String name = "EndToEnd";
		tempDir = Files.createTempDirectory(name);
		FileUtils.forceDeleteOnExit(tempDir.toFile());
		libraryFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/truncated_pan_human_library.dlib", tempDir, name, ".elib").toFile();
		libraryInterface = new LibraryFile() {{openFile(libraryFile);}};
		diaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia", tempDir, name, ".dia").toFile();

		diaFile2 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		diaFile3 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_3_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		fastaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tempDir, name, ".fasta").toFile();
		libraryScoringFactory = new EncyclopediaOneScoringFactory(parameters);
	}

	@After
	public void tearDown() throws Exception {
		parameters = null;
		if (null != libraryInterface) {
			libraryInterface.close();
			libraryInterface = null;
		}
		if (null != libraryFile){
			FileUtils.deleteQuietly(libraryFile);
			libraryFile = null;
		}
		if (null != diaFile) {
			FileUtils.deleteQuietly(diaFile);
			diaFile = null;
		}
		if (null != diaFile2) {
			FileUtils.deleteQuietly(diaFile2);
			diaFile2 = null;
		}
		if (null != diaFile3) {
			FileUtils.deleteQuietly(diaFile3);
			diaFile3 = null;
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
	public void testWholePipelineSingleData() throws Exception {
		EncyclopediaJobData jobDataA = new EncyclopediaJobData(diaFile,fastaFile,libraryInterface,libraryScoringFactory);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataA);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertSanityTest(outputFile,400,300);

		File tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();

		tempReport.delete();
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,400,300);
	}

	@Test
	public void testWholePipelineMultipleData() throws Exception {
		EncyclopediaJobData jobDataA = new EncyclopediaJobData(diaFile,fastaFile,libraryInterface,libraryScoringFactory);
		EncyclopediaJobData jobDataB = new EncyclopediaJobData(diaFile2,fastaFile,libraryInterface,libraryScoringFactory);
		EncyclopediaJobData jobDataC = new EncyclopediaJobData(diaFile3,fastaFile,libraryInterface,libraryScoringFactory);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataA);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertSanityTest(outputFile,400,300);

		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataB);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib")));

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib"));

		assertSanityTest(outputFile,400,300);

		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataC);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib")));

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib"));

		assertSanityTest(outputFile,400,400);

		File tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();

		tempReport.delete();
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,false);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,1200,400);
	}

	@Test
	public void testWholePipelineMultipleDataQuant() throws Exception {
		EncyclopediaJobData jobDataA = new EncyclopediaJobData(diaFile,fastaFile,libraryInterface,libraryScoringFactory);
		EncyclopediaJobData jobDataB = new EncyclopediaJobData(diaFile2,fastaFile,libraryInterface,libraryScoringFactory);
		EncyclopediaJobData jobDataC = new EncyclopediaJobData(diaFile3,fastaFile,libraryInterface,libraryScoringFactory);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataA);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataB);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataC);

		LibraryFile outputFile = new LibraryFile();

		//the output assertions here would be the same as for testWholePipelineMultipleData, no sense
		//checking them twice

		File tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();

		tempReport.delete();
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,1200,400);
	}

	public static void assertSanityTest(LibraryFile outputFile, int peptideFloor, int proteinFloor) throws Exception {
		assertTrue(MAX_POSSIBLE_PEPTIDES >= outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertTrue(peptideFloor <= outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertTrue(MAX_POSSIBLE_PROTEIN_GROUPS >= outputFile.getProteinGroups().size());
		assertTrue(proteinFloor <= outputFile.getProteinGroups().size());
		assertTrue(STANDARD_RANGE.contains(outputFile.getMinMaxMZ()));
	}
}
