package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.PecanParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
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

import static edu.washington.gs.maccoss.encyclopedia.EndToEndIT.assertSanityTest;
import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.assertTrue;

public class XCorDIAEndToEndIT {

	File diaFile;
	File diaFile2;
	File diaFile3;

	File fastaFile;

	File tempReport;

	PecanSearchParameters parameters;
	Path tempDir;

	XCorDIAOneScoringFactory factory;

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

		parameters = PecanParameterParser.getDefaultParametersObject();
		factory = new XCorDIAOneScoringFactory(parameters);
		String name = "EndToEnd";
		tempDir = Files.createTempDirectory(name);
		FileUtils.forceDeleteOnExit(tempDir.toFile());
		diaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia", tempDir, name, ".dia").toFile();

		diaFile2 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		diaFile3 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_3_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		fastaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tempDir, name, ".fasta").toFile();

		tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();
		tempReport.delete();
	}

	@After
	public void tearDown() throws Exception {
		parameters = null;
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
		if (null != tempReport){
			FileUtils.deleteQuietly(tempReport);
			tempReport = null;
		}
		if (null != tempDir) {
			FileUtils.deleteDirectory(tempDir.toFile());
			tempDir = null;
		}
	}

	@Test
	public void testWholePipelineSingleData() throws Exception {
		XCorDIAJobData jobDataA=new XCorDIAJobData(Optional.empty(), Optional.empty(), diaFile, fastaFile, new File(diaFile.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);

		XCorDIA.runPie(new EmptyProgressIndicator(), jobDataA);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".xcordia.txt")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertSanityTest(outputFile,400,300);
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,400,300);
	}

	@Test
	public void testWholePipelineMultipleData() throws Exception {
		XCorDIAJobData jobDataA=new XCorDIAJobData(Optional.empty(), Optional.empty(), diaFile, fastaFile, new File(diaFile.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);
		XCorDIAJobData jobDataB=new XCorDIAJobData(Optional.empty(), Optional.empty(), diaFile2, fastaFile, new File(diaFile2.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);
		XCorDIAJobData jobDataC=new XCorDIAJobData(Optional.empty(), Optional.empty(), diaFile3, fastaFile, new File(diaFile3.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);
		XCorDIA.runPie(new EmptyProgressIndicator(), jobDataA);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertSanityTest(outputFile,400,300);

		XCorDIA.runPie(new EmptyProgressIndicator(), jobDataB);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib")));

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib"));

		assertSanityTest(outputFile,400,300);

		XCorDIA.runPie(new EmptyProgressIndicator(), jobDataC);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib")));

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib"));

		assertSanityTest(outputFile,400,400);

		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,false);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,1200,400);
	}

	@Test
	public void testWholePipelineMultipleDataQuant() throws Exception {

		XCorDIAJobData jobDataA=new XCorDIAJobData(Optional.empty(), Optional.empty(), diaFile, fastaFile, new File(diaFile.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);
		XCorDIAJobData jobDataB=new XCorDIAJobData(Optional.empty(), Optional.empty(), diaFile2, fastaFile, new File(diaFile2.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);
		XCorDIAJobData jobDataC=new XCorDIAJobData(Optional.empty(), Optional.empty(), diaFile3, fastaFile, new File(diaFile3.getAbsolutePath()+XCorDIAJobData.OUTPUT_FILE_SUFFIX), factory);
		XCorDIA.runPie(new EmptyProgressIndicator(), jobDataA);
		XCorDIA.runPie(new EmptyProgressIndicator(),jobDataB);
		XCorDIA.runPie(new EmptyProgressIndicator(),jobDataC);

		LibraryFile outputFile = new LibraryFile();

		//the output assertions here would be the same as for testWholePipelineMultipleData, no sense
		//checking them twice
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,1200,400);
	}

}
