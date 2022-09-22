package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.base.Preconditions;
import com.google.common.collect.*;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.utils.ByteConverter;
import edu.washington.gs.maccoss.encyclopedia.utils.CompressionUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import junit.framework.AssertionFailedError;
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
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.List;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.*;

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
	static double UPPER_BOUND_PEPTIDE_MATCH = 1 / 0.85;
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

		try {
			libraryFile = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/truncated_pan_human_library.dlib", tempDir, name, ".dlib").toFile();
			libraryInterface = new LibraryFile();
			libraryInterface.openFile(libraryFile);

			diaFile = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia", tempDir, name, ".dia").toFile();
			diaFile2 = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia", tempDir, "EndToEnd", ".dia").toFile();
			diaFile3 = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_3_600.dia", tempDir, "EndToEnd", ".dia").toFile();

			fastaFile = getResourceAsTempFile(AbstractEndToEndIT.class, "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tempDir, name, ".fasta").toFile();
		} catch (AssumptionViolatedException exception) {
			LoggerFactory.getLogger(AbstractEndToEndIT.class)
					.warn("Assumption violated!", exception);
			throw exception;
		}
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
		SearchToBLIB.convert(new EmptyProgressIndicator(), ImmutableList.of(jobDataA,jobDataB,jobDataC), tempReport,false,true);
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

	/**
	 * Try running the two-step align-then-quantify workflow, quantifying only a subset of files.
	 * This test does not compare the results to a regression reference, but instead asserts that every peptide
	 * is quantified identically to the single-step align-and-quant approach used in {@link #testWholePipelineMultipleDataQuant()}.
	 */
	@Test
	public void testWholePipelineMultipleDataAlignOnlyQuantWorkflow() throws Exception {
		final ImmutableList<QuantitativeSearchJobData> jobData = ImmutableList.of(jobDataA, jobDataB, jobDataC);
		final SearchParameters parameters = jobData.iterator().next().getParameters();

		Assume.assumeTrue("Test requires quantitative search parameters", parameters.isQuantifySameFragmentsAcrossSamples());

		// First, run search with the normal single-step align-and-quant workflow as a reference comparison
		Logger.logLine("*********************************************");
		Logger.logLine("GENERATING REGULAR QUANT RESULTS AS REFERENCE");
		Logger.logLine("*********************************************");

		final Path standardQuantReport = Files.createTempFile(tempDir, "test_", ".elib");
		try {
			SearchToBLIB.convert(
					new EmptyProgressIndicator(),
					jobData,
					standardQuantReport.toFile(),
					SearchToBLIB.OutputFormat.ELIB,
					true,
					parameters
			);
		} catch (Exception e) {
			Assume.assumeNoException("Unable to generate reference results", e);
		}

		Logger.logLine("*********************************************");
		Logger.logLine("DONE WITH REGULAR QUANT RESULTS AS REFERENCE");
		Logger.logLine("*********************************************");

		Logger.logLine("*********************************************");
		Logger.logLine("CREATING ALIGNMENT-ONLY RESULTS FOR ALL JOBS");
		Logger.logLine("*********************************************");

		// Generate the alignment-only output
		SearchToBLIB.convert(new EmptyProgressIndicator(), jobData, tempReport, SearchToBLIB.OutputFormat.ALIB, true, parameters);

		Logger.logLine("*********************************************");
		Logger.logLine("DONE WITH ALIGNMENT-ONLY RESULTS FOR ALL JOBS");
		Logger.logLine("*********************************************");

		assertTrue("Output was created?", FileUtils.directoryContains(tempDir.toFile(), tempReport));

		// Perform some quick checks on the alignment-only results
		final LibraryFile outputFile = new LibraryFile();
		final int numAlignedPeptides;
		try {
			outputFile.openFile(tempReport);

			assertSanityTest(outputFile, getPeptideFloor(), getProteinFloor());

			try (Connection c = outputFile.getConnection()) {
				try (Statement s = c.createStatement()) {
					numAlignedPeptides = s.executeQuery("SELECT count() FROM entries").getInt(1);
				}
			}
		} finally {
			outputFile.close();
		}

		Logger.logLine("*********************************************");
		Logger.logLine("QUANTIFYING PEPTIDES IN A SUBSET OF JOBS");
		Logger.logLine("*********************************************");

		// Now we execute the quant step for two of the three jobs: the "seed" job and a non-seed, as RT alignment
		// differs between the two types of job.
		final ImmutableList<QuantitativeSearchJobData> quantJobData = jobData.subList(1, 3);
		final Path quantReport = Files.createTempFile(tempDir, "test_", ".elib");
		SearchToBLIB.convertElibQuantOnly(
				new EmptyProgressIndicator(),
				quantJobData,
				quantReport.toFile(),
				tempReport,
				quantJobData.iterator().next().getParameters()
		);

		Logger.logLine("*********************************************");
		Logger.logLine("DONE QUANTIFYING PEPTIDES IN A SUBSET OF JOBS");
		Logger.logLine("*********************************************");

		final LibraryFile quantFile = new LibraryFile();
		try {
			quantFile.openFile(quantReport.toFile());

			assertSanityTest(quantFile, getPeptideFloor(), 0); // not all might be quanted; no proteins are written

			// Check that both of the jobs have exactly the same results as run with the single-shot workflow.
			assertJobResultsMatch(quantJobData, quantFile, standardQuantReport);

			// This doesn't check just the quanted subset, so commented out for now.
			// Not that important because we've already made the stricter assertion that
			// this workflow matches the results of the default single-shot quant workflow.
//			assertValidBasedOnReference(quantFile, getReferenceMultiQuantResource());

			// Check that the non-quantified job(s) don't have results
			try (Connection c = quantFile.getConnection()) {
				try (PreparedStatement ps = c.prepareStatement(
						"SELECT count(), 'quants'" +
						" FROM peptidequants" +
						" WHERE sourcefile=?" +
						" UNION SELECT count(), 'rts'" +
						" FROM retentiontimes" +
						" WHERE sourcefile=?" +
						";"
				)) {
					final Set<SearchJobData> nonquantJobs = Sets.newHashSet(jobData);
					nonquantJobs.removeAll(quantJobData);

					for (SearchJobData job : nonquantJobs) {
						ps.setString(1, job.getDiaFileReader().getOriginalFileName());
						ps.setString(2, job.getDiaFileReader().getOriginalFileName());

						try (ResultSet rs = ps.executeQuery()) {
							assertTrue(rs.next());

							assertEquals("Too many " + rs.getString(2) + " for " + job.getDiaFileReader().getOriginalFileName() + ": " + rs.getInt(1),
									0,
									rs.getInt(1)
							);
						}
					}

				}
			}
		} finally {
			quantFile.close();
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

		assertTrue(
				String.format(
						"Fewer peptides than expected in %s: %d < %.02f",
						newFile.getName()+", (actual: "+peptides.size()+", expected:"+expectedPeptides.size()+")",
						peptides.size(),
						LOWER_BOUND_PEPTIDE_MATCH * expectedPeptides.size()
				),
				peptides.size() > LOWER_BOUND_PEPTIDE_MATCH * expectedPeptides.size()
		);

		assertTrue(
				String.format(
						"More peptides than expected in %s: %d > %.02f",
						newFile.getName()+", (actual: "+peptides.size()+", expected:"+expectedPeptides.size()+")",
						peptides.size(),
						UPPER_BOUND_PEPTIDE_MATCH * expectedPeptides.size()
				),
				peptides.size() < UPPER_BOUND_PEPTIDE_MATCH * expectedPeptides.size()
		);
		final long peptideMatches = peptides.stream()
				.filter(hasPeptideMatch(expectedPeptides))
				.count();

		// 85% of the peptides we IDed this run should be present in the previous results.
		// We don't bother checking if the same portion of the old results are still present,
		// this and the preceding checks for overall number are sufficiently reassuring.
		double percentage = peptideMatches / ((double) peptides.size());

		try {
			assertTrue("Fewer than 85% peptides match reference in " + newFile.getName()+", (actual: "+percentage+")", percentage > (1 / UPPER_BOUND_PEPTIDE_MATCH));
		} catch (AssertionError e) {
			// If the percentage is within epsilon of the bound, ignore the assertion failure.
			try {
				assertEquals(percentage, (1 / UPPER_BOUND_PEPTIDE_MATCH), 0.005f); // epsilon is 0.5%
			} catch (AssertionError e2) {
				throw e;
			}
		}

		assertTrue("pi0 lower than expected in " + newFile.getName()+", (actual: "+newFile.getMetadata().get("pi0")+")", Double.parseDouble(newFile.getMetadata().get("pi0")) > LOWER_BOUND_PI0_MATCH * (Double.parseDouble(expectedPi0)));
		assertTrue("pi0 greater than expected in " + newFile.getName()+", (actual: "+newFile.getMetadata().get("pi0")+")", Double.parseDouble(newFile.getMetadata().get("pi0")) < UPPER_BOUND_PI0_MATCH * (Double.parseDouble(expectedPi0)));
	}

	/**
	 * Check that the results for the given jobs match the given reference. More targeted version of {@link #assertValidBasedOnReference(LibraryFile, String)},
	 * but also much stricter -- expects the same set of quantifications, including identical RT and intensity.
	 *
	 * Also check associated quant reports exist and match.
	 */
	private void assertJobResultsMatch(Collection<? extends QuantitativeSearchJobData> jobs, LibraryFile quantFile, Path refElib) throws Exception {
		try (Connection c = quantFile.getConnection()) {
			try (PreparedStatement s = c.prepareStatement("ATTACH ? AS ref;")) {
				s.setString(1, refElib.toAbsolutePath().toString());

				s.execute();
			}

			final double epsilon = 1e-3;

			// Note: "seed" job will have no RT points, so this is commented
//			try (PreparedStatement s = c.prepareStatement(
//					"SELECT count()" +
//					" FROM retentiontimes" +
//					" WHERE SourceFile = ?;"
//			)) {
//				for (QuantitativeSearchJobData job : jobs) {
//					s.setString(1, job.getDiaFileReader().getOriginalFileName());
//
//					try (ResultSet rs = s.executeQuery()) {
//						assertTrue(rs.next());
//
//						final int numRts = rs.getInt(1);
//
//						Logger.logLine(String.format("Found %d RTs for %s in quant ELIB", numRts, job.getDiaFileReader().getOriginalFileName()));
//						assertNotEquals(0, numRts);
//					}
//				}
//			}

			// Check that each job's RT alignment data is the same (it should just be straight copied through)
			try (PreparedStatement s = c.prepareStatement(
					"SELECT t.peptidemodseq, t.predicted, r.predicted" +
					" FROM retentiontimes t" +
					" LEFT JOIN ref.retentiontimes r USING (peptidemodseq, sourcefile)" +
					" WHERE SourceFile = ?" +
					" AND (" +
							" abs(t.predicted - r.predicted) > ?" +
							" OR abs(t.predicted - r.predicted) > ?" +
					");"
			)) {
				for (QuantitativeSearchJobData job : jobs) {
					s.setString(1, job.getDiaFileReader().getOriginalFileName());
					s.setDouble(2, epsilon);
					s.setDouble(3, epsilon);

					try (ResultSet rs = s.executeQuery()) {
						while (rs.next()) {
							throw new AssertionError(String.format(
									"%s in %s: expected %.02f, got %.02f",
									rs.getString(1),
									job.getDiaFileReader().getOriginalFileName(),
									rs.getFloat(3),
									rs.getFloat(2)
							));
						}
					}
				}
			}

			// Check for quantified peptides not quantified in reference
			try (PreparedStatement s = c.prepareStatement(
					"SELECT count(), peptidemodseq" +
					" FROM peptidequants q" +
					" WHERE q.sourcefile = ?" +
					" AND NOT EXISTS (" +
					" SELECT 1 FROM ref.peptidequants rq" +
					" WHERE rq.peptidemodseq=q.peptidemodseq AND rq.precursorcharge=q.precursorcharge AND rq.sourcefile=q.sourcefile" +
					" LIMIT 1" +
					");"
			)) {
				for (QuantitativeSearchJobData job : jobs) {
					s.setString(1, job.getDiaFileReader().getOriginalFileName());

					try (ResultSet rs = s.executeQuery()) {
						assertTrue(rs.next());

						final int numMissingRef = rs.getInt(1);

						Logger.logLine(String.format("Found %d peptides not quantified in reference for %s, e.g. %s",
								numMissingRef,
								job.getDiaFileReader().getOriginalFileName(),
								rs.getString(2)
						));
						assertEquals(0, numMissingRef);
					}
				}
			}

			// Check for reference peptides not quantified
			try (PreparedStatement s = c.prepareStatement(
					"SELECT count(), q.peptidemodseq" +
					" FROM ref.peptidequants q" +
					" WHERE q.sourcefile = ?" +
					" AND NOT EXISTS (" +
					" SELECT 1 FROM peptidequants rq" +
					" WHERE rq.peptidemodseq=q.peptidemodseq AND rq.precursorcharge=q.precursorcharge AND rq.sourcefile=q.sourcefile" +
					" LIMIT 1" +
					");"
			)) {
				for (QuantitativeSearchJobData job : jobs) {
					s.setString(1, job.getDiaFileReader().getOriginalFileName());

					try (ResultSet rs = s.executeQuery()) {
						assertTrue(rs.next());

						final int numMissingQuant = rs.getInt(1);

						Logger.logLine(String.format("Found %d reference peptides not quantified for %s, e.g. %s",
								numMissingQuant,
								job.getDiaFileReader().getOriginalFileName(),
								rs.getString(2)
						));
						assertEquals(0, numMissingQuant);
					}
				}
			}

			// Check for different entries -- counts
			try (PreparedStatement s = c.prepareStatement(
					"SELECT count(), peptidemodseq" +
							" FROM entries e" +
							" LEFT JOIN ref.entries re USING (PeptideModSeq, PrecursorCharge, SourceFile)" +
							" WHERE SourceFile = ?" +
							" AND (" +
							" abs(e.rtinseconds - re.rtinseconds) > ?" +
							");"
			)) {
				for (QuantitativeSearchJobData job : jobs) {
					s.setString(1, job.getDiaFileReader().getOriginalFileName());
					s.setDouble(2, epsilon);

					try (ResultSet rs = s.executeQuery()) {
						assertTrue(rs.next());

						final int numMismatch = rs.getInt(1);

						Logger.logLine(String.format("Found %d mismatched entries for %s, e.g. %s",
								numMismatch,
								job.getDiaFileReader().getOriginalFileName(),
								rs.getString(2)
						));
//						assertEquals(0, numMismatch); // will be asserted row-by-row below
					}
				}
			}

			// Assert for each entry individually
			try (PreparedStatement s = c.prepareStatement(
					"SELECT" +
							" e.peptidemodseq," +
							" e.rtinseconds," +
							" re.rtinseconds" +
							" FROM entries e" +
							" LEFT JOIN ref.entries re USING (PeptideModSeq, PrecursorCharge, SourceFile)" +
							" WHERE SourceFile = ?;"
			)) {
				for (QuantitativeSearchJobData job : jobs) {
					s.setString(1, job.getDiaFileReader().getOriginalFileName());

					try (ResultSet rs = s.executeQuery()) {
						while (rs.next()) {
							final String pep = rs.getString(1);
							final float rt = rs.getFloat(2);
							final float refRt = rs.getFloat(3);

							assertEquals(
									String.format("entry rt mismatch: %s in %s: expected %.02f, got %.02f",
											pep,
											job.getDiaFileReader().getOriginalFileName(),
											refRt,
											rt
									),
									refRt,
									rt,
									epsilon
							);
						}
					}
				}
			}

			// Check for different quantification -- counts
			try (PreparedStatement s = c.prepareStatement(
					"SELECT count(), peptidemodseq" +
					" FROM peptidequants q" +
					" LEFT JOIN ref.peptidequants rq USING (PeptideModSeq, PrecursorCharge, SourceFile)" +
					" WHERE SourceFile = ?" +
					" AND (" +
							" abs(q.rtinsecondscenter - rq.rtinsecondscenter) > ?" +
							" OR abs(q.rtinsecondsstart - rq.rtinsecondsstart) > ?" +
							" OR abs(q.rtinsecondsstop - rq.rtinsecondsstop) > ?" +
							" OR abs(q.totalintensity - rq.totalintensity) > ?" +
					");"
			)) {
				for (QuantitativeSearchJobData job : jobs) {
					s.setString(1, job.getDiaFileReader().getOriginalFileName());
					s.setDouble(2, epsilon);
					s.setDouble(3, epsilon);
					s.setDouble(4, epsilon);
					s.setDouble(5, epsilon);

					try (ResultSet rs = s.executeQuery()) {
						assertTrue(rs.next());

						final int numMismatch = rs.getInt(1);

						Logger.logLine(String.format("Found %d mismatched quants for %s, e.g. %s",
								numMismatch,
								job.getDiaFileReader().getOriginalFileName(),
								rs.getString(2)
						));
//						assertEquals(0, numMismatch); // will be asserted row-by-row below
					}
				}
			}

			// Assert for each peptide's quant individually
			try (PreparedStatement s = c.prepareStatement(
					"SELECT" +
					" q.peptidemodseq," +
					" q.rtinsecondscenter," +
					" q.totalintensity," +
					" q.quantionmassarray," +
					" q.quantionmasslength," +
					" rq.rtinsecondscenter," +
					" rq.totalintensity," +
					" rq.quantionmassarray," +
					" rq.quantionmasslength" +
					" FROM peptidequants q" +
					" LEFT JOIN ref.peptidequants rq USING (PeptideModSeq, PrecursorCharge, SourceFile)" +
					" WHERE SourceFile = ?;"
			)) {
				for (QuantitativeSearchJobData job : jobs) {
					s.setString(1, job.getDiaFileReader().getOriginalFileName());

					try (ResultSet rs = s.executeQuery()) {
						while (rs.next()) {
							final String pep = rs.getString(1);
							final float rt = rs.getFloat(2);
							final float inten = rs.getFloat(3);
							final double[] ions = ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(4), rs.getInt(5)));
							final float refRt = rs.getFloat(6);
							final float refInten = rs.getFloat(7);
							final double[] refIons = ByteConverter.toDoubleArray(CompressionUtils.decompress(rs.getBytes(8), rs.getInt(9)));

							assertEquals(
									String.format("quant rt mismatch: %s in %s: expected (%.02f, %.02f), got (%.02f, %.02f)",
											pep,
											job.getDiaFileReader().getOriginalFileName(),
											refRt,
											refInten,
											rt,
											inten
									),
									refRt,
									rt,
									epsilon
							);
							assertArrayEquals(
									String.format("ions mismatch for %s in %s: expected %s, got %s",
											pep,
											job.getDiaFileReader().getOriginalFileName(),
											Arrays.toString(refIons),
											Arrays.toString(ions)
									),
									refIons,
									ions,
									epsilon
							);
							assertEquals(
									String.format("intensity mismatch for %s in %s: expected (%.02f, %.02f), got (%.02f, %.02f)",
											pep,
											job.getDiaFileReader().getOriginalFileName(),
											refRt,
											refInten,
											rt,
											inten
									),
									refInten,
									inten,
									epsilon
							);
						}
					}
				}
			}

			try (
					PreparedStatement s = c.prepareStatement(
							"SELECT p.peptideseq, p.proteinaccession" +
							" FROM peptidetoprotein p" +
							" LEFT JOIN ref.peptidetoprotein rp" +
							" USING (PeptideSeq, ProteinAccession)" +
							" WHERE rp.peptideseq IS NULL" +
							" AND NOT p.isdecoy" +
							" LIMIT 1;"
					);
					ResultSet rs = s.executeQuery()
			) {
				if (rs.next()) {
					throw new AssertionFailedError(String.format(
							"File has extra peptide-protein connections, e.g. (%s, %s)",
							rs.getString(1),
							rs.getString(2)
					));
				}
			}

			try (
					PreparedStatement s = c.prepareStatement(
							"SELECT rp.peptideseq, rp.proteinaccession" +
							" FROM ref.peptidetoprotein rp" +
							" LEFT JOIN peptidetoprotein p" +
							" USING (PeptideSeq, ProteinAccession)" +
							" WHERE p.peptideseq IS NULL" +
							" AND NOT p.isdecoy" +
							" LIMIT 1;"
					);
					ResultSet rs = s.executeQuery()
			) {
				if (rs.next()) {
					throw new AssertionFailedError(String.format(
							"File is missing some peptide-protein connections, e.g. (%s, %s)",
							rs.getString(1),
							rs.getString(2)
					));
				}
			}
		}

		SearchToBLIBIT.assertHasQuantReports(refElib);
		SearchToBLIBIT.assertHasQuantReports(quantFile.getFile().toPath());

		assertTableEquals(
				parseQuantTable(new File(refElib.toAbsolutePath() + ".peptides.txt"), "Peptide", jobs),
				parseQuantTable(new File(quantFile.getFile().getAbsolutePath() + ".peptides.txt"), "Peptide", jobs)
		);
		assertTableEquals(
				parseQuantTable(new File(refElib.toAbsolutePath() + ".proteins.txt"), "Protein", jobs),
				parseQuantTable(new File(quantFile.getFile().getAbsolutePath() + ".proteins.txt"), "Protein", jobs)
		);
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

	/**
	 * Parse the provided TSV quant report.
	 *
	 * @param rowKeyCol the column used as the key for each row
	 * @param jobData used to identify the columns that should be included
	 * @return rows: rowKey, cols: file name
	 */
	private static Table<String, String, Float> parseQuantTable(File f, String rowKeyCol, Collection<? extends SearchJobData> jobData) {
		final class QuantMuscle implements TableParserMuscle {
			private final Collection<String> quantCols;

			private final Table<String, String, Float> table = HashBasedTable.create();

			public QuantMuscle(Collection<String> quantCols) {
				this.quantCols = quantCols;
			}

			@Override
			public void processRow(Map<String, String> row) {
				final String rowKey = row.get(rowKeyCol);

				assertNotNull("No such row key in " + row.toString(), rowKey);

				for (String colKey : quantCols) {
					assertTrue("Missing column " + colKey + " in row " + rowKey, row.containsKey(colKey));

					table.put(rowKey, colKey, Float.parseFloat(row.get(colKey)));
				}
			}

			@Override
			public void cleanup() {
				// Nothing to do; must keep table in memory (for now)
			}

			Table<String, String, Float> finish() {
				final ImmutableTable<String, String, Float> result = ImmutableTable.copyOf(table);
				table.clear();
				return result;
			}
		}

		final QuantMuscle muscle = new QuantMuscle(jobData.stream().map(d -> d.getDiaFileReader().getOriginalFileName()).collect(Collectors.toSet()));

		TableParser.parseTSV(f, muscle);

		return muscle.finish();
	}

	/**
	 * More readable assertion of table equality; useful for large tables!
	 */
	private static void assertTableEquals(Table<String, String, Float> expected, Table<String, String, Float> actual) {
		expected.rowMap().forEach((rowKey, row) -> {
			final Map<String, Float> actualRow = actual.row(rowKey);

			assertNotNull("Did not find expected row with key " + rowKey, actualRow);
			assertFalse("Did not find expected row with key " + rowKey, actualRow.isEmpty());

			row.forEach((colKey, exp) -> {
				assertTrue(
						"Did not get intensity for " + rowKey + " in " + colKey,
						actualRow.containsKey(colKey)
				);
				assertEquals(
						"Intensity mismatches for " + rowKey + " in " + colKey,
						exp,
						actualRow.get(colKey),
						0.0001
				);
			});
		});

		// If there aren't the same number of rows the tables are unequal. If they do have the same size then the
		// check above ensures they're identical. Do this check last to give more useful feedback whenever possible.
		assertEquals("Wrong number of actual rows!", expected.rowKeySet().size(), actual.rowKeySet().size());

		assertFalse("No rows found in quant table!", actual.rowKeySet().isEmpty());
	}

	public static void assertSanityTest(LibraryFile outputFile, int peptideFloor, int proteinFloor) throws Exception {
		int peptideCount = outputFile.getAllEntries(false, AminoAcidConstants.createEmptyFixedAndVariable()).size();
		int proteinCount = outputFile.getProteinGroups().size();

		System.out.println("Peptides: " + peptideCount);
		System.out.println("Proteins: " + proteinCount);
		assertTrue("More than " + MAX_POSSIBLE_PEPTIDES + " peptides identified in " + outputFile.getName()+", (actual: "+peptideCount+")", MAX_POSSIBLE_PEPTIDES >= peptideCount);
		assertTrue("Fewer than " + peptideFloor + " peptides identified in " + outputFile.getName()+", (actual: "+peptideCount+")", peptideFloor <= peptideCount);
		assertTrue("More than " + MAX_POSSIBLE_PROTEIN_GROUPS + " protein groups identified in " + outputFile.getName()+", (actual: "+proteinCount+")", MAX_POSSIBLE_PROTEIN_GROUPS >= proteinCount);
		assertTrue("Fewer than " + proteinFloor + " protein groups identified in " + outputFile.getName()+", (actual: "+proteinCount+")", proteinFloor <= proteinCount);
		//assertTrue("Peptide identified outside of " + STANDARD_RANGE + " in " + outputFile.getName(), STANDARD_RANGE.contains(outputFile.getMinMaxMZ()));
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

		final Path resultsDir = targetDir.resolve("reference-data");

		final Path destination = Paths.get(resultsDir.toString(), resourcePath);
		Files.createDirectories(destination.getParent()); // ensure all folders under target/ exist

		Files.copy(
				elib.toPath(),
				destination,
				StandardCopyOption.REPLACE_EXISTING
		);
	}

	private static Path getElibFromResultsDirectory(String resourcePath) {
		return Paths.get("target", "reference-data", resourcePath);
	}
}
