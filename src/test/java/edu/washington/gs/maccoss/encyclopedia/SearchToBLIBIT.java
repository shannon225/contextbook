package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.QuantitativeSearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Assume;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.List;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;
import static org.junit.Assert.*;

public class SearchToBLIBIT {
	private static final String MOCK_PERCOLATOR_VERSION = "percolator_test_version";

	private final ProgressIndicator progress = new EmptyProgressIndicator();

	private SearchParameters searchParameters;
	private Path tempDir;

	private Path libraryA;
	private Path diaA;
	private Path featuresTxtA;
	private Path fastaA;
	private Path peptideOutputA;
	private Path decoyOutputA;

	@Before
	public void setUp() throws Exception {
		searchParameters = SearchParameterParser.getDefaultParametersObject();
		String name = "SearchToBLIBIT_";
		tempDir = Files.createTempDirectory(name);
		FileUtils.forceDeleteOnExit(tempDir.toFile());

		libraryA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/pan_human_library_600to603.dlib", tempDir, name, ".dlib");
		diaA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/bcs_2020jan16_600to603_hela_clib.dia", tempDir, name, ".dia");
		featuresTxtA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/bcs_2020jan16_600to603_hela_clib.dia.features.txt", tempDir, name, ".txt");
		fastaA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/pan_human_library_600to603.fasta", tempDir, name, ".fasta");
		peptideOutputA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/bcs_2020jan16_600to603_hela_clib.dia.encyclopedia.txt", tempDir, name, ".txt");
		decoyOutputA = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/testdata/bcs_2020jan16_600to603_hela_clib.dia.encyclopedia.decoy.txt", tempDir, name, ".txt");
	}

	@After
	public void tearDown() throws Exception {
		searchParameters = null;

		if (null != libraryA) {
			FileUtils.deleteQuietly(libraryA.toFile());
		}
		if (null != diaA) {
			FileUtils.deleteQuietly(diaA.toFile());
		}
		if (null != featuresTxtA) {
			FileUtils.deleteQuietly(featuresTxtA.toFile());
		}
		if (null != fastaA) {
			FileUtils.deleteQuietly(fastaA.toFile());
		}
		if (null != peptideOutputA) {
			FileUtils.deleteQuietly(peptideOutputA.toFile());
		}
		if (null != decoyOutputA) {
			FileUtils.deleteQuietly(decoyOutputA.toFile());
		}
		if (null != tempDir) {
			FileUtils.deleteDirectory(tempDir.toFile());
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

	private SearchJobData getSearchJobDataA() throws IOException, SQLException {
		return makeJobData(libraryA, diaA, featuresTxtA, fastaA, peptideOutputA, decoyOutputA);
	}

	private SearchJobData getSearchJobDataB() throws IOException, SQLException {
		Assume.assumeTrue(false); //TODO: get data B!!
		throw new UnsupportedOperationException("TODO");
//		//TODO: move to resources
//		final Path dia = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/LongSwath_UPS1_40fm_Ecoli_1ug-rep2.dia");
//		final Path featuresTxt = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep2.dia.features.txt");
//		final Path fasta = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/ups-protein-standards.fasta");
//		final Path peptideOutput = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep2.dia.xcordia.txt");
//		final Path decoyOutput = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep2.dia.xcordia.decoy.txt");
//
//		return makeXCorDIAJobData(dia, featuresTxt, fasta, peptideOutput, decoyOutput);
	}

	private QuantitativeSearchJobData makeJobData(Path library, Path dia, Path featuresTxt, Path fasta, Path peptideOutput, Path decoyOutput) throws IOException, SQLException {
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
		);
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