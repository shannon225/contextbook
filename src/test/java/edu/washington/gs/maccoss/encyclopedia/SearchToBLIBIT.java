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

import java.io.File;
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
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertTrue("Result file had no entries", 0 < numEntries);

		assertHasPercolatorMetadata(file);

		//TODO: other assertions specific to this output format
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