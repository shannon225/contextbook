package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.junit.AfterClass;
import org.junit.AssumptionViolatedException;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.sql.Connection;
import java.sql.Statement;
import java.util.HashMap;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class EncyclopediaEndToEndIT extends AbstractEndToEndIT{
	static final String REFERENCE_SEARCH1_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/encyc_dia_1.elib";
	static final String REFERENCE_SEARCH2_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/encyc_dia_2.elib";
	static final String REFERENCE_SEARCH3_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/encyc_dia_3.elib";
	static final String REFERENCE_SINGLE_QUANT_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/encyc_single.elib";
	static final String REFERENCE_MULTI_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/encyc_multi.elib";
	static final String REFERENCE_MULTI_QUANT_RESOURCE = "/edu/washington/gs/maccoss/encyclopedia/reference_data/encyc_multi_quant.elib";

	static SearchParameters parameters;
	static LibraryScoringFactory libraryScoringFactory;

	static int PEPTIDE_FLOOR = 400;
	static int PROTEIN_FLOOR = 300;

	@BeforeClass
	public static void buildReports() throws Exception {
		// IMPORTANT: we must run the search with the same parameters as the reference data.
		// This can be challenging when the default parameters change! While we can hard-code
		// certain settings, it can be hard to figure out the changes and copy them here,
		// so instead we fetch the reference parameters directly. Note that a side effect of this
		// will be that new versions of reference artifacts built with this code will use whatever
		// parameters the configured reference used.

		final HashMap<String, String> metadata;

		final LibraryFile reference = new LibraryFile();
		try {
			final Path refElib = EncyclopediaTestUtils.getResourceAsTempFile(EncyclopediaEndToEndIT.class, REFERENCE_SEARCH1_RESOURCE, tempDir, "reference_", ".elib");
			reference.openFile(refElib.toFile());
			metadata = Maps.newHashMap(reference.getMetadata());
		} finally {
			reference.close();
		}

		// This must be true for most tests, regardless of what the reference did.
		metadata.put("-quantifyAcrossSamples", "true");

		parameters = SearchParameterParser.parseParameters(metadata);

		libraryScoringFactory = new EncyclopediaOneScoringFactory(parameters);
		setUpClass();
		jobDataA = makeAndDoJob(diaFile);
		jobDataB = makeAndDoJob(diaFile2);
		jobDataC = makeAndDoJob(diaFile3);

		copyJobDataToResultsDirectory(jobDataA, REFERENCE_SEARCH1_RESOURCE);
		copyJobDataToResultsDirectory(jobDataB, REFERENCE_SEARCH2_RESOURCE);
		copyJobDataToResultsDirectory(jobDataC, REFERENCE_SEARCH3_RESOURCE);
	}

	@AfterClass
	public static void tearDownReports() throws Exception {
		tearDownClass();
		parameters = null;
		libraryScoringFactory = null;
	}

	public static EncyclopediaJobData makeAndDoJob(File dia) throws Exception {
		EncyclopediaJobData jobData = new EncyclopediaJobData(dia,fastaFile,libraryInterface,libraryScoringFactory);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobData);
		return jobData;
	}

	@Override
	public int getPeptideFloor() {
		return PEPTIDE_FLOOR;
	}

	@Override
	public int getProteinFloor() {
		return PROTEIN_FLOOR;
	}

	@Override
	public List<String> getReferenceSearchResources() throws Exception {
		return ImmutableList.of(
				REFERENCE_SEARCH1_RESOURCE,
				REFERENCE_SEARCH2_RESOURCE,
				REFERENCE_SEARCH3_RESOURCE
		);
	}

	@Override
	public String getReferenceSingleQuantResource() throws Exception {
		return REFERENCE_SINGLE_QUANT_RESOURCE;
	}

	@Override
	public String getReferenceMultiResource() throws Exception {
		return REFERENCE_MULTI_RESOURCE;
	}

	@Override
	public String getReferenceMultiQuantResource() throws Exception {
		return REFERENCE_MULTI_QUANT_RESOURCE;
	}

	/**
	 * Test running a search using the "alignment only" library (ALIB) output by the `-alignOnly` mode.
	 * This is <emph>not a supported mode of operation</emph>, but may happen accidentally. It will likely not give
	 * exactly the same results as searching a more normally-generated library, but should succeed and
	 * at least give reasonable results.
	 */
	@Test
	public void testWholePipelineMultipleDataQuantUsingAlignmentOnlyLibrary() throws Exception {
		final ImmutableList<QuantitativeSearchJobData> initialJobData = ImmutableList.of(jobDataA, jobDataB, jobDataC);
		final List<SearchJobData> jobData = Lists.newArrayList();

		final SearchParameters parameters = libraryScoringFactory.getParameters();

		final LibraryFile alignmentLib = new LibraryFile();
		final int numAlignedPeptides;
		try {
			// Generate the alignment-only output
			final Path alignmentFile = Files.createTempFile(tempDir, "test_", ".elib");
			try {
				SearchToBLIB.convert(new EmptyProgressIndicator(), initialJobData, alignmentFile.toFile(), SearchToBLIB.OutputFormat.ALIB, true, parameters);

				assertTrue("Output was created?", Files.exists(alignmentFile));

				alignmentLib.openFile(alignmentFile.toFile());

				// Perform some quick checks on the alignment-only results
				assertSanityTest(alignmentLib, getPeptideFloor(), getProteinFloor());

				try (Connection c = alignmentLib.getConnection()) {
					try (Statement s = c.createStatement()) {
						numAlignedPeptides = s.executeQuery("SELECT count() FROM entries").getInt(1);
					}
				}
			} catch (Exception e) {
				throw new AssumptionViolatedException("Unable to generate reference results", e);
			}

			// Now we run the "normal" search-and-quant approach, but using the alignment-only library.
			// This requires regenerating the job data.
			for (QuantitativeSearchJobData job : initialJobData) {
				final Path newDia = Files.createTempFile(tempDir, "test_", ".dia");
				Files.copy(job.getDiaFileReader().getFile().toPath(), newDia, StandardCopyOption.REPLACE_EXISTING);

				final EncyclopediaJobData newJob = new EncyclopediaJobData(
						newDia.toFile(),
						fastaFile,
						alignmentLib,
						libraryScoringFactory
				);
				jobData.add(newJob);

				Encyclopedia.runSearch(new EmptyProgressIndicator(), newJob);
			}

			SearchToBLIB.convert(
					new EmptyProgressIndicator(),
					jobData,
					tempReport,
					SearchToBLIB.OutputFormat.ELIB,
					true,
					parameters
			);
		} finally {
			alignmentLib.close();
		}

		final LibraryFile quantFile = new LibraryFile();
		try {
			quantFile.openFile(tempReport);

			// Only require 2/3rds of the ID count of the normal search (this is a degraded library)
			assertSanityTest(quantFile, (int) Math.round(0.67 * getPeptideFloor() * jobData.size()), (int) Math.round(0.67 * getProteinFloor()));

			// We don't assert any validity based on reference, as long as we get sane numbers we consider this
			// to be "good enough".
//			assertValidBasedOnReference(quantFile, getReferenceMultiQuantResource());
		} finally {
			quantFile.close();
		}
	}
}
