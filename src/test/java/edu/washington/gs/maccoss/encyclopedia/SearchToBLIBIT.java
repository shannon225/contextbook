package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.AlternatePeakLocationInferrer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutor;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableConcatenator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Assume;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.DataFormatException;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.*;

public class SearchToBLIBIT {
	private static final String MOCK_PERCOLATOR_VERSION = "percolator_test_version";
	static final double DELTA = 0.0001;

	private final ProgressIndicator progress = new EmptyProgressIndicator();

	private SearchParameters searchParameters;
	private Path tempDir;

	private Path library;
	private Path fasta;

	private Path diaA;
	private Path elibA;
	private Path featuresTxtA;
	private Path peptideOutputA;
	private Path decoyOutputA;

	private Path diaB;
	private Path elibB;
	private Path featuresTxtB;
	private Path peptideOutputB;
	private Path decoyOutputB;

	@Before
	public void setUp() throws Exception {
		searchParameters = SearchParameterParser.getDefaultParametersObject();
		String name = "SearchToBLIBIT_";
		tempDir = Files.createTempDirectory(name);
		FileUtils.forceDeleteOnExit(tempDir.toFile());

		library = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/truncated_pan_human_library.dlib", tempDir, name, ".dlib");
		fasta = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/uniprot_human_2018.subset.fasta", tempDir, name, ".fasta");

		diaA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia", tempDir, name, ".dia");
		elibA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia.elib", tempDir, name, ".elib");
		featuresTxtA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia.features.txt", tempDir, name, ".txt");
		peptideOutputA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia.encyclopedia.txt", tempDir, name, ".txt");
		decoyOutputA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_1_600.dia.encyclopedia.decoy.txt", tempDir, name, ".txt");

		diaB = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia", tempDir, name, ".dia");
		elibB = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia.elib", tempDir, name, ".elib");
		featuresTxtB = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia.features.txt", tempDir, name, ".txt");
		peptideOutputB = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia.encyclopedia.txt", tempDir, name, ".txt");
		decoyOutputB = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/121115_bcs_hela_24mz_400_1000_0D_2_600.dia.encyclopedia.decoy.txt", tempDir, name, ".txt");
	}

	@After
	public void tearDown() throws Exception {
		searchParameters = null;

		if (null != library) {
			FileUtils.deleteQuietly(library.toFile());
			library = null;
		}
		if (null != fasta) {
			FileUtils.deleteQuietly(fasta.toFile());
			fasta = null;
		}

		if (null != diaA) {
			FileUtils.deleteQuietly(diaA.toFile());
			diaA = null;
		}
		if (null != featuresTxtA) {
			FileUtils.deleteQuietly(featuresTxtA.toFile());
			featuresTxtA = null;
		}
		if (null != peptideOutputA) {
			FileUtils.deleteQuietly(peptideOutputA.toFile());
			peptideOutputA = null;
		}
		if (null != decoyOutputA) {
			FileUtils.deleteQuietly(decoyOutputA.toFile());
			decoyOutputA = null;
		}

		if (null != diaB) {
			FileUtils.deleteQuietly(diaB.toFile());
			diaB = null;
		}
		if (null != featuresTxtB) {
			FileUtils.deleteQuietly(featuresTxtB.toFile());
			featuresTxtB = null;
		}
		if (null != peptideOutputB) {
			FileUtils.deleteQuietly(peptideOutputB.toFile());
			peptideOutputB = null;
		}
		if (null != decoyOutputB) {
			FileUtils.deleteQuietly(decoyOutputB.toFile());
			decoyOutputB = null;
		}

		if (null != tempDir) {
			FileUtils.deleteDirectory(tempDir.toFile());
			tempDir = null;
		}
	}

	@Test
	public void testConvertMultiSampleElib() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(
				getSearchJobDataA(),
				getSearchJobDataB()
		);

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false, // elib
				true
		);

		final LibraryFile file = new LibraryFile();
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertTrue("Result file had no entries", 0 < numEntries);

		assertHasPercolatorMetadata(file);
	}

	@Test
	public void testConvertMultiSampleBlib() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".blib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(
				getSearchJobDataA(),
				getSearchJobDataB()
		);

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				true, // blib
				true
		);

		assertValidBlib(libFile);
	}

	@Test
	public void testConvertMultiSampleQuant() throws Exception {
		// create quant parameters
		final HashMap<String, String> parameterMap = searchParameters.toParameterMap();
		parameterMap.put("-quantifyAcrossSamples", "true");
		searchParameters = SearchParameterParser.parseParameters(parameterMap);

		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(
				getSearchJobDataA(),
				getSearchJobDataB()
		);

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false, // elib
				true
		);

		final LibraryFile file = new LibraryFile();
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertTrue("Result file had no entries", 0 < numEntries);

		assertHasPercolatorMetadata(file);
	}

	@Test
	public void testConvertMultiSampleAlignOnly() throws Exception {
		// create quant parameters
		final HashMap<String, String> parameterMap = searchParameters.toParameterMap();
		parameterMap.put("-quantifyAcrossSamples", "true");
		searchParameters = SearchParameterParser.parseParameters(parameterMap);

		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(
				getSearchJobDataA(),
				getSearchJobDataB()
		);

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				SearchToBLIB.OutputFormat.ALIB,
				true
		);

		final LibraryFile file = new LibraryFile();
		try {
			file.openFile(libFile.toFile());

			final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
			assertTrue("Result file had no entries", 0 < numEntries);

			assertValidAlib(file, jobData);
		} finally {
			file.close();
		}
	}

	/**
	 * Same as above, but at a lower level to enable direct comparison of the inferrer.
	 */
	@Test
	public void testConvertMultiSampleAlignOnlyAlignment() throws Exception {
		// create quant parameters
		final HashMap<String, String> parameterMap = searchParameters.toParameterMap();
		parameterMap.put("-quantifyAcrossSamples", "true");
		searchParameters = SearchParameterParser.parseParameters(parameterMap);

		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(
				getSearchJobDataA(),
				getSearchJobDataB()
		);

		// We reproduce the behavior of this method to get access to lower-level info, i.e. `inferrer`
//		SearchToBLIB.convert(progress,
//				jobData,
//				libFile.toFile(),
//				SearchToBLIB.OutputFormat.ALIB,
//				true
//		);
		final SearchJobData representativeJob = jobData.iterator().next();
		String filename=libFile.toFile().getName();
		if (filename.lastIndexOf('.')>0) {
			filename=filename.substring(0, filename.lastIndexOf('.'));
		}
		File bigFeatureFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_features.txt");
		File bigPercolatorFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_results.txt");
		File bigPercolatorDecoyFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_decoy.txt");
		File bigPercolatorProteinFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_protein_results.txt");
		File bigPercolatorProteinDecoyFile=new File(representativeJob.getPercolatorFiles().getInputTSV().getParentFile(), filename+"_concatenated_protein_decoy.txt");
		PercolatorExecutionData bigPercolatorFiles=new PercolatorExecutionData(bigFeatureFile, representativeJob.getPercolatorFiles().getFastaFile(), bigPercolatorFile, bigPercolatorDecoyFile, bigPercolatorProteinFile, bigPercolatorProteinDecoyFile, searchParameters);

		Logger.logLine("Running global Percolator analysis.");
		TableConcatenator.concatenateTables(
				jobData.stream()
						.map(SearchJobData::getPercolatorFiles)
						.map(PercolatorExecutionData::getInputTSV)
						.collect(Collectors.toCollection(ArrayList::new)),
				bigFeatureFile
		);

		int modelNumber = Integer.MAX_VALUE; // always use the last model (if reusing a model)
		final Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides = PercolatorExecutor.executePercolatorTSV(searchParameters.getPercolatorVersionNumber(), bigPercolatorFiles, searchParameters.getEffectivePercolatorThreshold(), searchParameters.getAAConstants(), modelNumber);

		final PeakLocationInferrerInterface inferrer = AlternatePeakLocationInferrer.getAlignmentData(new EmptyProgressIndicator(), jobData, passingPeptides.x, searchParameters);

		SearchToBLIB.OutputFormat.ALIB.convert(
				new EmptyProgressIndicator(),
				jobData,
				libFile.toFile(),
				passingPeptides,
				Optional.of(bigPercolatorFiles),
				Optional.of(inferrer),
				searchParameters
		);

		final LibraryFile file = new LibraryFile();
		try {
			file.openFile(libFile.toFile());

			final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
			assertTrue("Result file had no entries", 0 < numEntries);

			assertValidAlib(file, jobData);

			assertSameInferrer(inferrer, SearchToBLIB.readInferrer(file, jobData, searchParameters), jobData);
		} finally {
			file.close();
		}
	}

	/**
	 * Test what happens running the quant-only conversion ({@code -alignmentFrom}).
	 * This is essentially an extension of the method above, as we generate the alignment-only
	 * file, then separately quantify a single file.
	 */
	@Test
	public void testConvertMultiSampleQuantOnly() throws Exception {
		// create quant parameters
		final HashMap<String, String> parameterMap = searchParameters.toParameterMap();
		parameterMap.put("-quantifyAcrossSamples", "true");
		searchParameters = SearchParameterParser.parseParameters(parameterMap);

		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final Path quantFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(quantFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(quantFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(
				getSearchJobDataA(),
				getSearchJobDataB()
		);

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				SearchToBLIB.OutputFormat.ALIB,
				true
		);

		final LibraryFile file = new LibraryFile();
		try {
			file.openFile(libFile.toFile());

			final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
			assertTrue("Result file had no entries", 0 < numEntries);

			assertValidAlib(file, jobData);
		} catch (AssertionError e) {
			Assume.assumeNoException("Test setup failed: unable to produce valid alignment-only results", e);
		} finally {
			file.close();
		}

		SearchToBLIB.convertElibQuantOnly(progress,
//				jobData.subList(0, 1), // only the first job -- SUCCEEDS (this is the "seed" file)
				jobData.subList(1, 2), // only the second job -- more difficult b/c RT alignment is involved
				quantFile.toFile(),
				libFile.toFile(),
				searchParameters
		);

		final LibraryFile quantLib = new LibraryFile();
		try {
			quantLib.openFile(quantFile.toFile());

			try (Connection c = quantLib.getConnection()) {
				try (Statement s = c.createStatement()) {
					try (ResultSet rs = s.executeQuery("SELECT count() FROM peptidequants;")) {
						assertTrue(rs.next());
						assertTrue("Result file had no quants", rs.getInt(1) > 0);
					}
				}
			}

			//TODO: other assertions specific to this output format
		} finally {
			quantLib.close();
		}
	}

	/**
	 * Test what happens running the quant-only conversion ({@code -alignmentFrom})
	 * with a "normal" search library rather than an alignment-only results file.
	 */
	@Test(expected = NullPointerException.class) // failure due to ELIB missing pi0
	public void testConvertMultiSampleQuantOnlyWithNormalElib() throws Exception {
		// create quant parameters
		final HashMap<String, String> parameterMap = searchParameters.toParameterMap();
		parameterMap.put("-quantifyAcrossSamples", "true");
		searchParameters = SearchParameterParser.parseParameters(parameterMap);

		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(
				getSearchJobDataA()
		);

		SearchToBLIB.convertElibQuantOnly(progress,
				jobData,
				libFile.toFile(),
				// Read alignment from the input (search) library
				((LibraryFile) ((EncyclopediaJobData) jobData.iterator().next()).getLibrary()).getFile(),
				searchParameters
		);

		final LibraryFile file = new LibraryFile();
		try {
			file.openFile(libFile.toFile());

			try (Connection c = file.getConnection()) {
				try (Statement s = c.createStatement()) {
					try (ResultSet rs = s.executeQuery("SELECT count() FROM peptidequants;")) {
						assertTrue(rs.next());
						assertTrue("Result file had no quants", rs.getInt(1) > 0);
					}
				}
			}

			//TODO: other assertions specific to this output format
		} finally {
			file.close();
		}
	}

	@Test
	public void testConvertSingleSampleElib() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(getSearchJobDataA());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false, // elib
				true
		);

		final LibraryFile file = new LibraryFile();
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertTrue("Result file had no entries", 0 < numEntries);

		assertHasPercolatorMetadata(file);
		assertEquals("Found unexpected Percolator version in output ELIB", MOCK_PERCOLATOR_VERSION, file.getMetadata().get(LibraryFile.PERCOLATOR_VERSION));
	}

	@Test
	public void testConvertSingleSampleBlib() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".blib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(getSearchJobDataA());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				true, // blib
				true
		);

		assertValidBlib(libFile);
	}

	@Test
	public void testConvertSingleSampleQuant() throws Exception {
		// create quant parameters
		final HashMap<String, String> parameterMap = searchParameters.toParameterMap();
		parameterMap.put("-quantifyAcrossSamples", "true");
		searchParameters = SearchParameterParser.parseParameters(parameterMap);

		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = ImmutableList.of(getSearchJobDataA());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false, // elib
				true
		);

		final LibraryFile file = new LibraryFile();
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertTrue("Result file had no entries", 0 < numEntries);

		assertHasPercolatorMetadata(file);
		assertEquals("Found unexpected Percolator version in output ELIB", MOCK_PERCOLATOR_VERSION, file.getMetadata().get(LibraryFile.PERCOLATOR_VERSION));
	}

	private void assertValidBlib(Path blib) throws IOException {
		assertTrue("BLIB doesn't exist!", Files.exists(blib));

		assertTrue("BLIB is too short!", 1024L < Files.size(blib));
	}

	private void assertHasPercolatorMetadata(LibraryFile file) throws IOException, SQLException {
		if (null == file) {
			throw new NullPointerException("Can't run assertions on null LibraryFile!");
		}

		final HashMap<String, String> metadata = file.getMetadata();

		metadata.forEach((k, v) -> System.out.println(String.format("%s:\t%s", k, v)));

		assertNotNull(metadata.get(LibraryFile.PERCOLATOR_VERSION));
		assertNotNull(metadata.get("pi0"));
	}

	/**
	 * Check that the file appears to be a valid and correctly-written ALIB.
	 *
	 * @param jobData jobs to check for in the output
	 *
	 * @see SearchToBLIB.OutputFormat#ALIB
	 */
	private void assertValidAlib(LibraryFile file, List<SearchJobData> jobData) throws SQLException, IOException, DataFormatException {
		try (Connection c = file.getConnection()) {
			try (Statement s = c.createStatement()) {
				try (ResultSet rs = s.executeQuery("SELECT count() FROM entries;")) {
					assertTrue(rs.next());
					assertTrue("ALIB had no peptide entries!", 10 < rs.getInt(1));
				}

				try (ResultSet rs = s.executeQuery("SELECT count() FROM peptidescores;")) {
					assertTrue(rs.next());
					assertTrue("ALIB had no scored peptides!", 10 < rs.getInt(1));
				}

				try (ResultSet rs = s.executeQuery("SELECT count() FROM peptidetoprotein;")) {
					assertTrue(rs.next());
					assertTrue("ALIB had no peptide-protein connections!", 10 < rs.getInt(1));
				}

				final Pair<ArrayList<PercolatorPeptide>, Float> passingPeptides = SearchToBLIB.readPassingPeptides(file, searchParameters);
				assertTrue("Unable to fetch passing peptides from ALIB!", 0 < passingPeptides.x.size());

				try (ResultSet rs = s.executeQuery("SELECT PeptideModSeq, count(), massencodedlength/8 FROM entries GROUP BY PeptideModSeq;")) {
					while (rs.next()) {
						final String pep = rs.getString(1);
						final int nEntries = rs.getInt(2);
						final int nIons = rs.getInt(3);

						assertEquals(pep + " had wrong number of entries: " + nEntries, 1, nEntries);

						// Note that these entries aren't filtered  for the min number of quant ions;
						// that filter happens during quantification.
						assertTrue(pep + " had unexpected number of quant fragments: " + nIons,
								nIons >= 0 || nIons <= searchParameters.getNumberOfQuantitativePeaks()
						);
					}
				}
			}

			try (PreparedStatement ps = c.prepareStatement(
					"SELECT count() FROM retentiontimes WHERE sourcefile = ?;"
			)) {
				boolean foundSeed = false;
				for (SearchJobData job : jobData) {
					ps.setString(1, job.getDiaFileReader().getOriginalFileName());

					try (ResultSet rs = ps.executeQuery()) {
						assertTrue(rs.next());

						final int nRtPoints = rs.getInt(1);

						if (nRtPoints < 10) {
							assertFalse("Already found seed job! Not enough RT points (" + nRtPoints + ") for " + job.getDiaFileReader().getOriginalFileName(), foundSeed);
							foundSeed = true;
						}
					}
				}
			}

			try (PreparedStatement ps = c.prepareStatement(
					"SELECT count()" +
					" FROM retentiontimes rt" +
					" JOIN entries e USING (peptidemodseq)" +
					" WHERE abs(e.rtinseconds - rt.library) > 0.001;"
			)) {
				try (ResultSet rs = ps.executeQuery()) {
					assertTrue(rs.next());

					final int nEntries = rs.getInt(1);

					assertEquals(
							"Too many mismatched RTs! All entries should match recorded \"library\" RT!",
							0,
							nEntries
					);
				}
			}

			// Same check as above, but with the decoded `inferrer`.
			final PeakLocationInferrerInterface inferrer = SearchToBLIB.readInferrer(file, jobData, searchParameters);
			for (SearchJobData job : jobData) {
				final List<RetentionTimeAlignmentInterface.AlignmentDataPoint> alignmentData = inferrer.getAlignmentData(job);
				// Just ensure the inferred RTs match the corresponding saved points.
				boolean hadPeptidePoint = false;
				for (RetentionTimeAlignmentInterface.AlignmentDataPoint p : alignmentData) {
					if (null == p.getPeptideModSeq()) {
						// In an ALIB not all RT points will be for a peptide, some may just be alignment "knots"
						continue;
					} else if (!hadPeptidePoint) {
						hadPeptidePoint = true;
					}

					assertEquals(
							String.format("%s in %s",
									p.getPeptideModSeq(),
									job.getDiaFileReader().getOriginalFileName()
							),
							p.getPredictedActual(), // in mins
							inferrer.getWarpedRTInSec(job, p.getPeptideModSeq()) / 60f, // convert to mins
							DELTA
					);
				}

				if (!hadPeptidePoint) {
					// This is the "seed" job, so we need different assertions -- check that inferrer's
					// warped RT matches the saved alignment RT (derived from this job).
					for (LibraryEntry e : file.getAllEntries(false, searchParameters.getAAConstants())) {
						assertEquals(
								String.format("%s in %s",
										e.getPeptideModSeq(),
										job.getDiaFileReader().getOriginalFileName()
								),
								e.getRetentionTime(), // in sec
								inferrer.getWarpedRTInSec(job, e.getPeptideModSeq()),
								DELTA
						);
					}
				}
			}
		}

		assertHasPercolatorMetadata(file);
	}

	private void assertSameInferrer(PeakLocationInferrerInterface expected, PeakLocationInferrerInterface actual, List<SearchJobData> jobs) {
		final Set<String> allPeptides = jobs.stream()
				.flatMap(j -> expected.getAlignmentData(j).stream())
				.map(RetentionTimeAlignmentInterface.AlignmentDataPoint::getPeptideModSeq)
				.collect(Collectors.toSet());

		// Check each peptide's predicted RT in each job
		for (SearchJobData job : jobs) {
			for (String modSeq : allPeptides) {
				assertEquals(
						String.format("wrong warped RT: %s in %s",
								modSeq,
								job.getDiaFileReader().getOriginalFileName()
						),
						expected.getWarpedRTInSec(job, modSeq),
						actual.getWarpedRTInSec(job, modSeq),
						DELTA
				);
			}
		}
	}

	private SearchJobData getSearchJobDataA() throws IOException, SQLException {
		return makeJobData(library, diaA, featuresTxtA, fasta, peptideOutputA, decoyOutputA, elibA);
	}

	private SearchJobData getSearchJobDataB() throws IOException, SQLException {
		return makeJobData(library, diaB, featuresTxtB, fasta, peptideOutputB, decoyOutputB, elibB);
	}

	private QuantitativeSearchJobData makeJobData(Path library, Path dia, Path featuresTxt, Path fasta, Path peptideOutput, Path decoyOutput, Path resultsElib) throws IOException, SQLException {
		Assume.assumeTrue(Files.exists(dia));

		final StripeFile diaReader = new StripeFile(true) ;
		diaReader.openFile(dia.toFile());

		final TestPercolatorExecutionData percolatorFiles = new TestPercolatorExecutionData(featuresTxt, fasta, peptideOutput, decoyOutput);

		// Set up the state as though we've just generated these files using Percolator.
		percolatorFiles.setPercolatorExecutableVersion(MOCK_PERCOLATOR_VERSION);

		return new EncyclopediaJobData(
				dia.toFile(), // dia file; must exist
				diaReader,
				percolatorFiles,
				searchParameters,
				"TEST",
				new LibraryFile() {{ openFile(library.toFile()); }},
				new EncyclopediaOneScoringFactory(searchParameters)
		) {
			@Override
			public File getResultLibrary() {
				// Must ensure we grab the right temp file instead of using path munging
				return resultsElib.toFile();
			}
		};
	}

	/**
	 * This subclass allows access to the {@link #setPercolatorExecutableVersion(String)}
	 */
	private class TestPercolatorExecutionData extends PercolatorExecutionData {
		public TestPercolatorExecutionData(Path featuresTxt, Path fasta, Path peptideOutput, Path decoyOutput) {
			super(
					featuresTxt.toFile(),
					fasta.toFile(),
					peptideOutput.toFile(),
					decoyOutput.toFile(),
					null,
					null,
					SearchToBLIBIT.this.searchParameters
			);
		}

		@Override
		public void setPercolatorExecutableVersion(String percolatorExecutableVersion) {
			super.setPercolatorExecutableVersion(percolatorExecutableVersion);
		}

	}
}