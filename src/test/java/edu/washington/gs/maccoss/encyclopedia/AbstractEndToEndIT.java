package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultimap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.junit.*;
import org.slf4j.LoggerFactory;
import java.util.Collection;
import java.util.List;

import java.awt.*;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import java.util.function.Predicate;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.assertTrue;

public abstract class AbstractEndToEndIT {

	static File diaFile;
	static File diaFile2;
	static File diaFile3;

	static File fastaFile;
	static File libraryFile;
	static LibraryFile libraryInterface;

	static File tempReport;

	static Path tempDir;

	static Range STANDARD_RANGE = new Range(592.5840338877389,604.3740813086648);
	static int MAX_POSSIBLE_PEPTIDES = 4669;
	static int MAX_POSSIBLE_PROTEIN_GROUPS = 6676;

	@BeforeClass
	public static void buildReports() throws Exception {
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
		libraryFile = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/truncated_pan_human_library.dlib", tempDir, name, ".dlib").toFile();
		libraryInterface = new LibraryFile() {{openFile(libraryFile);}};
		diaFile = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia", tempDir, name, ".dia").toFile();

		diaFile2 = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		diaFile3 = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_3_600.dia", tempDir, "EndToEnd", ".dia").toFile();
		fastaFile = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tempDir, name, ".fasta").toFile();

		tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();
		tempReport.delete();
	}

	@Before
	public void setUp() throws Exception {

	}

	@After
	public void tearDown() throws Exception {
		if (null != tempReport){
			FileUtils.deleteQuietly(tempReport);
		}
	}

	@AfterClass
	public static void tearDownReports() throws Exception {
		tempReport = null;
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
	}

	@Test
	public void testWholePipelineSingleData() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib")));

		LibraryFile outputFile = new LibraryFile();
		outputFile.openFile(FileUtils.getFile(tempDir.toFile(),diaFile.getName() + ".elib"));

		assertSanityTest(outputFile,getPeptideFloor(),getProteinFloor());
		assertValidBasedOnReference(outputFile,getReferenceSearches()[0]);
	}

	@Test
	public void testWholePipelineSingleDataQuant() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		LibraryFile outputFile = new LibraryFile();
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,getPeptideFloor(),getProteinFloor());

		assertValidBasedOnReference(outputFile,getReferenceSingleQuant());
	}

	@Test
	public void testWholePipelineMultipleData() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		SearchJobData jobDataB = makeAndDoJob(diaFile2);
		SearchJobData jobDataC = makeAndDoJob(diaFile3);

		LibraryFile outputFile = new LibraryFile();
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,false);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,getPeptideFloor() * 3,getProteinFloor());

		assertValidBasedOnReference(outputFile,getReferenceMulti());
	}

	@Test
	public void testWholePipelineMultipleDataQuant() throws Exception {
		SearchJobData jobDataA = makeAndDoJob(diaFile);
		SearchJobData jobDataB = makeAndDoJob(diaFile2);
		SearchJobData jobDataC = makeAndDoJob(diaFile3);

		LibraryFile outputFile = new LibraryFile();
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		outputFile.openFile(tempReport);

		assertSanityTest(outputFile,getPeptideFloor() * 3,getProteinFloor());

		assertValidBasedOnReference(outputFile,getReferenceMultiQuant());
	}

	public static void assertValidBasedOnReference(LibraryFile newFile, LibraryFile reference) throws Exception {
		List<LibraryEntry> peptides = newFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable());
		List<LibraryEntry> expectedPeptides = reference.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable());

		assertTrue (peptides.size() > (0.95) * expectedPeptides.size()
				&& peptides.size() < (1.05) * expectedPeptides.size());

		final long peptideMatches = peptides.stream()
				.filter(hasPeptideMatch(expectedPeptides))
				.count();

		// 95% of the peptides we IDed this run should be present in the previous results.
		// We don't bother checking if 95% of the old results are still present, this and
		// the precedingchecks for overall number are sufficiently reassuring.
		double percentage = peptideMatches / ((double) peptides.size());

		assertTrue(percentage > 0.95);

		assertTrue(Double.parseDouble(newFile.getMetadata().get("pi0")) > (0.75) * (Double.parseDouble(reference.getMetadata().get("pi0")))
				&& Double.parseDouble(newFile.getMetadata().get("pi0")) < (1.25) * (Double.parseDouble(reference.getMetadata().get("pi0"))));
	}

	/**
	 * @return A predicate that returns true if an element of {@code expectedPeptides} has
	 *         the same source and {@code peptideModSeq}, and satisfies {@link #isRtMatch(LibraryEntry)}.
	 */
	private static Predicate<? super LibraryEntry> hasPeptideMatch(Collection<LibraryEntry> expectedPeptides) {
		final ImmutableMultimap.Builder<String, LibraryEntry> b = ImmutableMultimap.builder();
		expectedPeptides.forEach(e -> b.put(e.getPeptideModSeq(), e));

		final ImmutableMultimap<String, LibraryEntry> expectedPeptidesByModSeq = b.build();

		return entry -> expectedPeptidesByModSeq.get(entry.getPeptideModSeq()).stream()
				.anyMatch(
						isRtMatch(entry)
						.and(e -> entry.getSource().equals(e.getSource())) // must be in same file
				);
	}

	/**
	 * @return A predicate that returns true if both the given entry and {@code entry}
	 *         are both {@link ChromatogramLibraryEntry} instances and have sufficiently
	 *         overlapping RT ranges.
	 */
	private static Predicate<LibraryEntry> isRtMatch(LibraryEntry entry) {
		Preconditions.checkArgument(entry instanceof ChromatogramLibraryEntry);

		final Range rtRange = ((ChromatogramLibraryEntry) entry).getRtRange();

		return e2 -> {
			Preconditions.checkState(e2 instanceof ChromatogramLibraryEntry);

			final Range r2 = ((ChromatogramLibraryEntry) e2).getRtRange();

			//TODO: assess degree of overlap
			return r2.getStart() < rtRange.getStop()
					&& r2.getStop() > rtRange.getStart();
		};
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
