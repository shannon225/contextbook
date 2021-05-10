package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultimap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.junit.*;
import org.slf4j.LoggerFactory;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.assertTrue;

public abstract class AbstractEndToEndIT {
	/**
	 * When making regression checks, what proportion of the reference result's
	 * RT range must be covered by the current result?
	 */
	static final double REQUIRED_PROPORTION_OF_REFERENCE_RT_RANGE = 0.60;

	static File diaFile;
	static File diaFile2;
	static File diaFile3;

	static File fastaFile;
	static File libraryFile;
	static LibraryFile libraryInterface;

	static File tempReport;

	static Path tempDir;

	static Range STANDARD_RANGE = new Range(592.5840338877389,604.3740813086648);

	static double LOWER_BOUND_PEPTIDE_MATCH = 0.95;
	static double UPPER_BOUND_PEPTIDE_MATCH = 1.05;
	static double LOWER_BOUND_PI0_MATCH = 0.75;
	static double UPPER_BOUND_PI0_MATCH = 1.25;

	static QuantitativeSearchJobData jobDataA;
	static QuantitativeSearchJobData jobDataB;
	static QuantitativeSearchJobData jobDataC;

	/**
	 * These values are the total peptides and proteing groups defined by
	 * the libraryFile used in the test.
	 */
	static int MAX_POSSIBLE_PEPTIDES = 4669;
	static int MAX_POSSIBLE_PROTEIN_GROUPS = 6676;

	@BeforeClass
	public static void setUpClass() throws Exception {
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
	}

	@Before
	public void setUp() throws Exception {
		tempReport = Files.createTempFile(tempDir, "test_",".elib").toFile();
		tempReport.delete();
	}

	@After
	public void tearDown() throws Exception {
		if (null != tempReport){
			FileUtils.deleteQuietly(tempReport);
			tempReport = null;
		}
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
		tempReport = null;
		jobDataA = null;
		jobDataB = null;
		jobDataC = null;
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
		final File resultLibrary = jobDataA.getResultLibrary();

		assertTrue(FileUtils.directoryContains(tempDir.toFile(), resultLibrary));

		final LibraryFile outputFile = new LibraryFile();
		try {
			outputFile.openFile(resultLibrary);

			final String referenceResource = getReferenceSearchResources().get(0);

			// This job's ELIB will already have been copied to the target directory during setup

			assertSanityTest(outputFile, getPeptideFloor(), getProteinFloor());
			assertValidBasedOnReference(outputFile, referenceResource);
		} finally {
			outputFile.close();
		}
	}

	@Test
	public void testWholePipelineSingleDataQuant() throws Exception {
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA),tempReport,false,true);

		assertTrue(FileUtils.directoryContains(tempDir.toFile(), tempReport));

		LibraryFile outputFile = new LibraryFile();
		try {
			outputFile.openFile(tempReport);

			final String referenceResource = getReferenceSingleQuantResource();

			// Copy the data before checking assertions
			copyElibToResultsDirectory(tempReport, referenceResource);

			assertSanityTest(outputFile, getPeptideFloor(), getProteinFloor());
			assertValidBasedOnReference(outputFile, referenceResource);
		} finally {
			outputFile.close();
		}
	}

	@Test
	public void testWholePipelineMultipleData() throws Exception {
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,false);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		LibraryFile outputFile = new LibraryFile();
		try {
			outputFile.openFile(tempReport);

			final String referenceResource = getReferenceMultiResource();

			// Copy the data before checking assertions
			copyElibToResultsDirectory(tempReport, referenceResource);

			assertSanityTest(outputFile, getPeptideFloor() * 3, getProteinFloor());
			assertValidBasedOnReference(outputFile, referenceResource);
		} finally {
			outputFile.close();
		}
	}

	@Test
	public void testWholePipelineMultipleDataQuant() throws Exception {
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC),tempReport,false,true);
		assertTrue(FileUtils.directoryContains(tempDir.toFile(),tempReport));

		LibraryFile outputFile = new LibraryFile();
		try {
			outputFile.openFile(tempReport);

			final String referenceResource = getReferenceMultiQuantResource();

			// Copy the data before checking assertions
			copyElibToResultsDirectory(tempReport, referenceResource);

			assertSanityTest(outputFile, getPeptideFloor() * 3, getProteinFloor());
			assertValidBasedOnReference(outputFile, referenceResource);
		} finally {
			outputFile.close();
		}
	}

	public static void assertValidBasedOnReference(LibraryFile newFile, String referenceResource) throws Exception {
		final List<LibraryEntry> expectedPeptides;
		final String expectedPi0;
		{
			final Path refElib = getResourceAsTempFile(AbstractEndToEndIT.class, referenceResource, tempDir, "reference_", ".elib");
			final LibraryFile f = new LibraryFile();
			try {
				f.openFile(refElib.toFile());

				expectedPeptides = f.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable());
				expectedPi0 = f.getMetadata().get("pi0");
			} finally {
				f.close();
				FileUtils.deleteQuietly(refElib.toFile());
			}
		}

		final List<LibraryEntry> peptides = newFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable());

		assertTrue ("Fewer peptides than expected in " + newFile.getName(), peptides.size() > LOWER_BOUND_PEPTIDE_MATCH * expectedPeptides.size());
		assertTrue("More peptides than expected in " + newFile.getName(), peptides.size() < UPPER_BOUND_PEPTIDE_MATCH * expectedPeptides.size());

		final long peptideMatches = peptides.stream()
				.filter(hasPeptideMatch(expectedPeptides))
				.count();

		// 95% of the peptides we IDed this run should be present in the previous results.
		// We don't bother checking if 95% of the old results are still present, this and
		// the preceding checks for overall number are sufficiently reassuring.
		double percentage = peptideMatches / ((double) peptides.size());

		assertTrue("Fewer than 95% peptides match reference in " + newFile.getName(), percentage > LOWER_BOUND_PEPTIDE_MATCH);

		assertTrue("pi0 lower than expected in " + newFile.getName(), Double.parseDouble(newFile.getMetadata().get("pi0")) > LOWER_BOUND_PI0_MATCH * (Double.parseDouble(expectedPi0)));
	    assertTrue("pi0 greater than expected in " + newFile.getName(), Double.parseDouble(newFile.getMetadata().get("pi0")) < UPPER_BOUND_PI0_MATCH * (Double.parseDouble(expectedPi0)));
	}

	/**
	 * @return A predicate that returns true if an element of {@code expectedPeptides} has
	 *         the same source and {@code peptideModSeq}, and satisfies {@link #isRtMatch(LibraryEntry)}.
	 */
	private static Predicate<? super LibraryEntry> hasPeptideMatch(Collection<LibraryEntry> expectedPeptides) {
		final ImmutableMultimap.Builder<String, LibraryEntry> b = ImmutableMultimap.builder();
		expectedPeptides.forEach(e -> b.put(e.getPeptideModSeq(), e));

		final ImmutableMultimap<String, LibraryEntry> expectedPeptidesByModSeq = b.build();

		return entry -> {
			if (!(entry instanceof ChromatogramLibraryEntry)){
				//TODO: Figure out why non-ChromatogramLibraryEntry items appear in the library
				return false;
			}
			return expectedPeptidesByModSeq.get(entry.getPeptideModSeq()).stream()
					.filter( p -> p instanceof ChromatogramLibraryEntry)
					.anyMatch(
							isRtMatch(entry)
									.and(e -> entry.getSource().equals(e.getSource())) // must be in same file
					);
		};
	}

	/**
	 * @return A predicate that returns true if both the given entry and {@code entry}
	 *         are both {@link ChromatogramLibraryEntry} instances and have sufficiently
	 *         overlapping RT ranges.
	 */
	private static Predicate<LibraryEntry> isRtMatch(LibraryEntry entry) {
		Preconditions.checkArgument(entry instanceof ChromatogramLibraryEntry);

		final com.google.common.collect.Range<Float> guavaRange;
		{
			final Range rtRange = ((ChromatogramLibraryEntry) entry).getRtRange();
			guavaRange = com.google.common.collect.Range.closed(rtRange.getStart(), rtRange.getStop());
		}

		return reference -> {
			Preconditions.checkState(reference instanceof ChromatogramLibraryEntry);

			final Range r2 = ((ChromatogramLibraryEntry) reference).getRtRange();

			final com.google.common.collect.Range<Float> intersection;
			try {
				intersection = guavaRange.intersection(com.google.common.collect.Range.closed(r2.getStart(), r2.getStop()));
			} catch (IllegalArgumentException e){
				return false;
			}

			return
					intersection.hasLowerBound()
					&& intersection.hasUpperBound()

					// check that the intersection is at least 60% of the _reference entry_'s width
					&& intersection.upperEndpoint() - intersection.lowerEndpoint() > REQUIRED_PROPORTION_OF_REFERENCE_RT_RANGE * r2.getRange()
			;
		};
	}

	public static void assertSanityTest(LibraryFile outputFile, int peptideFloor, int proteinFloor) throws Exception {
		int peptideCount = outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size();
		int proteinCount = outputFile.getProteinGroups().size();

		System.out.println("Peptides: " + peptideCount);
		System.out.println("Proteins: " + proteinCount);
		assertTrue("More than " + MAX_POSSIBLE_PEPTIDES + " peptides identified in " + outputFile.getName(), MAX_POSSIBLE_PEPTIDES >= peptideCount);
		assertTrue("Fewer than " + peptideFloor + " peptides identified in " + outputFile.getName(), peptideFloor <= peptideCount);
		assertTrue("More than " + MAX_POSSIBLE_PROTEIN_GROUPS + " protein groups identified in " + outputFile.getName(), MAX_POSSIBLE_PROTEIN_GROUPS >= proteinCount);
		assertTrue("Fewer than " + proteinFloor + " protein groups identified in " + outputFile.getName(), proteinFloor <= proteinCount);
		assertTrue("Peptide identified outside of " + STANDARD_RANGE + " in " + outputFile.getName(), STANDARD_RANGE.contains(outputFile.getMinMaxMZ()));
	}

	public abstract int getPeptideFloor();

	public abstract int getProteinFloor();

	public abstract List<String> getReferenceSearchResources() throws Exception;

	public abstract String getReferenceSingleQuantResource() throws Exception;

	public abstract String getReferenceMultiResource() throws Exception;

	public abstract String getReferenceMultiQuantResource() throws Exception;

	protected static void copyJobDataToResultsDirectory(SearchJobData jobData, String resourcePath) throws IOException {
		copyElibToResultsDirectory(((QuantitativeSearchJobData) jobData).getResultLibrary(), resourcePath);
	}

	private static void copyElibToResultsDirectory(File elib, String resourcePath) throws IOException {
		final Path targetDir = Paths.get("target");
		Assume.assumeTrue(Files.exists(targetDir));

		final Path resultsDir = targetDir.resolve("regression-test-data");

		final Path destination = Paths.get(resultsDir.toString(), resourcePath);
		Files.createDirectories(destination.getParent()); // ensure all folders under target/ exist

		Files.copy(
				elib.toPath(),
				destination,
				StandardCopyOption.REPLACE_EXISTING
		);
	}
}
