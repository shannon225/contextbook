package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.gui.framework.SearchToELIBJob;
import edu.washington.gs.maccoss.encyclopedia.gui.general.JobProcessorTableModel;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.ForkJoinPool;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.*;

public class EndToEndIT {

	File diaFile;
	File fastaFile;
	LibraryInterface libraryInterface;
	LibraryScoringFactory libraryScoringFactory;
	SearchParameters parameters;
	Path tempDir;

	JobProcessorTableModel tableModel;

	@Before
	public void setUp() throws Exception {
		parameters = SearchParameterParser.getDefaultParametersObject();
		String name = "EndToEnd";
		tempDir = Files.createTempDirectory(name);
		FileUtils.forceDeleteOnExit(tempDir.toFile());

		libraryInterface = BlibToLibraryConverter.getFile(getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/truncated_pan_human_library.dlib", tempDir, name, ".elib").toFile());
		diaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia", tempDir, name, ".dia").toFile();
		fastaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tempDir, name, ".fasta").toFile();
		libraryScoringFactory = new EncyclopediaOneScoringFactory(parameters);

		tableModel = new JobProcessorTableModel();
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
		tableModel = null;
	}

	@Test
	public void testWholePipelineSingleData() throws Exception {
		EncyclopediaJobData jobDataA = new EncyclopediaJobData(diaFile,fastaFile,libraryInterface,libraryScoringFactory);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataA);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertEquals(407,outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertEquals(348,outputFile.getProteinGroups().size());
		assertEquals(new Range(592.6138073613809,604.3740813086648),outputFile.getMinMaxMZ());

		File tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();

		SearchToELIBJob job=new SearchToELIBJob(tempReport, true, tableModel);
		tableModel.addJob(job);

		tempReport.delete();
		assertNotNull(job);
		SearchToBLIB.convert(job.getProgressIndicator(), ImmutableList.of(jobDataA),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertEquals(407,outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertEquals(348,outputFile.getProteinGroups().size());
		assertEquals(new Range(592.6138073613809,604.3740813086648),outputFile.getMinMaxMZ());
	}

	@Test
	public void testWholePipelineMultipleData() throws Exception {
		File diaFile2 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		File diaFile3 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_3_600.dia", tempDir, "EndToEnd", ".dia").toFile();

		EncyclopediaJobData jobDataA = new EncyclopediaJobData(diaFile,fastaFile,libraryInterface,libraryScoringFactory);
		EncyclopediaJobData jobDataB = new EncyclopediaJobData(diaFile2,fastaFile,libraryInterface,libraryScoringFactory);
		EncyclopediaJobData jobDataC = new EncyclopediaJobData(diaFile3,fastaFile,libraryInterface,libraryScoringFactory);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataA);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertEquals(407,outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertEquals(348,outputFile.getProteinGroups().size());
		assertEquals(new Range(592.6138073613809,604.3740813086648),outputFile.getMinMaxMZ());

		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataB);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib")));

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib"));

		assertEquals(450,outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertEquals(374,outputFile.getProteinGroups().size());
		assertEquals(new Range(592.5840338877389,604.3740813086648),outputFile.getMinMaxMZ());

		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobDataC);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib")));

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib"));

		assertEquals(554,outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertEquals(458,outputFile.getProteinGroups().size());
		assertEquals(new Range(592.6138073613809,604.3740813086648),outputFile.getMinMaxMZ());

		File tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();

		SearchToELIBJob job=new SearchToELIBJob(tempReport, false, tableModel);
		tableModel.addJob(job);

		tempReport.delete();
		assertNotNull(job);
		SearchToBLIB.convert(job.getProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,false);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertEquals(1890,outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertEquals(481,outputFile.getProteinGroups().size());
		assertEquals(new Range(592.5840338877389,604.3740813086648),outputFile.getMinMaxMZ());
	}

	@Test
	public void testWholePipelineMultipleDataQuant() throws Exception {
		File diaFile2 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		File diaFile3 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_3_600.dia", tempDir, "EndToEnd", ".dia").toFile();

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

		SearchToELIBJob job=new SearchToELIBJob(tempReport, true, tableModel);
		tableModel.addJob(job);

		tempReport.delete();
		assertNotNull(job);
		SearchToBLIB.convert(job.getProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertEquals(1716,outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size());
		assertEquals(481,outputFile.getProteinGroups().size());
		assertEquals(new Range(592.5840338877389,604.3740813086648),outputFile.getMinMaxMZ());
	}
}
