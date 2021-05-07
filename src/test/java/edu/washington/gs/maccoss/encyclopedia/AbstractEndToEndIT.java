package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.LoggerFactory;
import java.util.List;

import java.awt.*;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.assertTrue;

public abstract class AbstractEndToEndIT {

	File diaFile;
	File diaFile2;
	File diaFile3;

	File fastaFile;
	File libraryFile;
	LibraryFile libraryInterface;

	File tempReport;

	Path tempDir;

	static Range STANDARD_RANGE = new Range(592.5840338877389,604.3740813086648);
	static int MAX_POSSIBLE_PEPTIDES = 4669;
	static int MAX_POSSIBLE_PROTEIN_GROUPS = 6676;

	@Before
	public void setUp() throws Exception {
		if (GraphicsEnvironment.isHeadless() && "1.8".equals(System.getProperty("java.specification.version"))) {
			LoggerFactory.getLogger(EncyclopediaEndToEndIT.class)
					.info("Disabling assistive technologies to avoid errors in headless build on Java 8!");

			// On JDK 8 running headless we can encounter problems if this is set by the system
			// installation; by overwriting it we avoid errors with JFreeChart. This is not an
			// issue on later releases that will better respect the headless flag.
			// See https://stackoverflow.com/a/59397731/115714
			System.setProperty("javax.accessibility.assistive_technologies", "java.lang.Object");
		}

		String name = "EndToEnd";
		tempDir = Files.createTempDirectory(name);
		FileUtils.forceDeleteOnExit(tempDir.toFile());
		libraryFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/truncated_pan_human_library.dlib", tempDir, name, ".dlib").toFile();
		libraryInterface = new LibraryFile() {{openFile(libraryFile);}};
		diaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia", tempDir, name, ".dia").toFile();

		diaFile2 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		diaFile3 = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_3_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		fastaFile = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tempDir, name, ".fasta").toFile();

		tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();
		tempReport.delete();
	}

	@After
	public void tearDown() throws Exception {
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
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertSanityTest(outputFile,getPeptideFloor(),getProteinFloor());
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,getPeptideFloor(),getProteinFloor());
	}

	@Test
	public void testWholePipelineMultipleData() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		SearchJobData jobDataB = makeAndDoJob(diaFile2);
		SearchJobData jobDataC = makeAndDoJob(diaFile3);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertSanityTest(outputFile,getPeptideFloor(),getProteinFloor());

		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib")));

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib"));

		assertSanityTest(outputFile,getPeptideFloor(),getProteinFloor());

		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib")));

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib"));

		assertSanityTest(outputFile,getPeptideFloor(),getProteinFloor());

		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,false);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,getPeptideFloor() * 3,getProteinFloor());
	}

	@Test
	public void testWholePipelineMultipleDataQuant() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		SearchJobData jobDataB = makeAndDoJob(diaFile2);
		SearchJobData jobDataC = makeAndDoJob(diaFile3);

		LibraryFile outputFile = new LibraryFile();

		//the output assertions here would be the same as for testWholePipelineMultipleData, no sense
		//checking them twice
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,getPeptideFloor() * 3,getProteinFloor());
	}

	@Test
	public void testSingleDataRegression() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertValidBasedOnReference(outputFile,getReferenceSearches()[0]);
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertValidBasedOnReference(outputFile,getReferenceSingleQuant());
	}

	@Test
	public void testMultipleDataRegression() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		SearchJobData jobDataB = makeAndDoJob(diaFile2);
		SearchJobData jobDataC = makeAndDoJob(diaFile3);

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));
		assertValidBasedOnReference(outputFile,getReferenceSearches()[0]);

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile2.getName() + ".elib"));
		assertValidBasedOnReference(outputFile,getReferenceSearches()[1]);

		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile3.getName() + ".elib"));
		assertValidBasedOnReference(outputFile,getReferenceSearches()[2]);

		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,false);
		outputFile.openFile(tempReport);
		assertValidBasedOnReference(outputFile,getReferenceMulti());
	}

	@Test
	public void testMultipleDataQuantRegression() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		SearchJobData jobDataB = makeAndDoJob(diaFile2);
		SearchJobData jobDataC = makeAndDoJob(diaFile3);

		LibraryFile outputFile = new LibraryFile();
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,false);
		outputFile.openFile(tempReport);
		assertValidBasedOnReference(outputFile,getReferenceMultiQuant());
	}

	public static void assertValidBasedOnReference(LibraryFile newFile, LibraryFile reference) throws Exception {
		List<LibraryEntry> peptides = newFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable());
		List<LibraryEntry> expectedPeptides = reference.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable());

		System.out.println(peptides.size());
		System.out.println(expectedPeptides.size());

		assertTrue (peptides.size() > (0.95) * expectedPeptides.size()
				&& peptides.size() < (1.05) * expectedPeptides.size());

		int peptideMatches = 0;

		for (LibraryEntry entry : peptides) {
			if (expectedPeptides.stream().anyMatch(e ->
					e.getPeptideModSeq().equals(entry.getPeptideModSeq())
					&& (Math.abs(e.getRetentionTime() - entry.getRetentionTime()) < 0.01))){
				peptideMatches++;
			}
		}

		double percentage = peptideMatches / ((double)Math.min(peptides.size(),expectedPeptides.size()));

		assertTrue(percentage > 0.95);
	}

	public static void assertSanityTest(LibraryFile outputFile, int peptideFloor, int proteinFloor) throws Exception {
		int peptideCount = outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size();
		int proteinCount = outputFile.getProteinGroups().size();

		System.out.println("Peptides: " + peptideCount);
		System.out.println("Proteins: " + proteinCount);
		assertTrue(MAX_POSSIBLE_PEPTIDES >= peptideCount);
		assertTrue(peptideFloor <= peptideCount);
		assertTrue(MAX_POSSIBLE_PROTEIN_GROUPS >= proteinCount);
		assertTrue(proteinFloor <= proteinCount);
		assertTrue(STANDARD_RANGE.contains(outputFile.getMinMaxMZ()));
	}

	public abstract SearchJobData makeAndDoJob(File dia) throws Exception;

	public abstract int getPeptideFloor();

	public abstract int getProteinFloor();

	public abstract LibraryFile[] getReferenceSearches() throws Exception;

	public abstract LibraryFile getReferenceSingleQuant() throws Exception;

	public abstract LibraryFile getReferenceMulti() throws Exception;

	public abstract LibraryFile getReferenceMultiQuant() throws Exception;
}
