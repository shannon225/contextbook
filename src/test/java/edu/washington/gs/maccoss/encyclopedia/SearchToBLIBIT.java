package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import org.apache.commons.io.FileUtils;
import org.apache.tools.ant.taskdefs.optional.extension.LibFileSet;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.mockito.internal.matchers.Null;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class SearchToBLIBIT {
	private static final int NUM_ENTRIES = 8;

	private final ProgressIndicator progress = new EmptyProgressIndicator();

	private boolean previousOpenInPlace;
	private SearchParameters searchParameters;
	private Path tempDir;

	@Before
	public void setUp() throws Exception {
		previousOpenInPlace = LibraryFile.OPEN_IN_PLACE;
		LibraryFile.OPEN_IN_PLACE = true; // non-default

		searchParameters = SearchParameterParser.getDefaultParametersObject();

		tempDir = Files.createTempDirectory("SearchToBLIBIT_");
		FileUtils.forceDeleteOnExit(tempDir.toFile());
	}

	@After
	public void tearDown() throws Exception {
		LibraryFile.OPEN_IN_PLACE = previousOpenInPlace; // restore default

		searchParameters = null;

		if (null != tempDir) {
			FileUtils.deleteDirectory(tempDir.toFile());
		}
	}

	@Test
	public void testConvertMultiSampleElib() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = Stream.of("test1", "test2").parallel()
				.map(this::createMockJobData)
				.collect(Collectors.toList());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false, // elib
				true
		);

		final LibraryFile file = new LibraryFile();
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertEquals(NUM_ENTRIES, numEntries); //TODO

		assertHasPercolatorMetadata(file);
	}

	@Test
	public void testConvertMultiSampleBlib() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".blib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = Stream.of("test1", "test2").parallel()
				.map(this::createMockJobData)
				.collect(Collectors.toList());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				true, // blib
				true
		);

		//TODO: assertions for blib
//		final LibraryFile file = new LibraryFile();
//		file.openFile(libFile.toFile());
//
//		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
//		assertEquals(NUM_ENTRIES, numEntries);
//
//		assertHasPercolatorMetadata(file);
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

		final List<SearchJobData> jobData = Stream.of("test1", "test2").parallel()
				.map(this::createMockJobData)
				.collect(Collectors.toList());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false, // elib
				true
		);

		final LibraryFile file = new LibraryFile();
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertEquals(NUM_ENTRIES, numEntries); //TODO

		assertHasPercolatorMetadata(file);
	}

	@Test
	public void testConvertSingleSampleElib() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".elib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = Stream.of("test1")
				.map(this::createMockJobData)
				.collect(Collectors.toList());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false, // elib
				true
		);

		final LibraryFile file = new LibraryFile();
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertEquals(NUM_ENTRIES, numEntries); //TODO

		assertHasPercolatorMetadata(file);
	}

	@Test
	public void testConvertSingleSampleBlib() throws Exception {
		final Path libFile = Files.createTempFile(tempDir, "SearchToBLIBIT_", ".blib");
		Files.delete(libFile); // can't exist (we're trying to create it)
		FileUtils.forceDeleteOnExit(libFile.toFile());

		final List<SearchJobData> jobData = Stream.of("test1")
				.map(this::createMockJobData)
				.collect(Collectors.toList());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				true, // blib
				true
		);

		//TODO: assertions for blib
//		final LibraryFile file = new LibraryFile();
//		file.openFile(libFile.toFile());
//
//		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
//		assertEquals(NUM_ENTRIES, numEntries);
//
//		assertHasPercolatorMetadata(file);
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

		final List<SearchJobData> jobData = Stream.of("test1")
				.map(this::createMockJobData)
				.collect(Collectors.toList());

		SearchToBLIB.convert(progress,
				jobData,
				libFile.toFile(),
				false, // elib
				true
		);

		final LibraryFile file = new LibraryFile();
		file.openFile(libFile.toFile());

		final int numEntries = file.getAllEntries(false, searchParameters.getAAConstants()).size();
		assertEquals(NUM_ENTRIES, numEntries); //TODO

		assertHasPercolatorMetadata(file);
	}

	private EncyclopediaJobData createMockJobData(String name) {
		try {
			final Path peptideOutput = Files.createTempFile(tempDir, name, ".peptides.txt");
			try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(peptideOutput.toFile())))) {
				pw.println("PSMId\tq-value\tposterior_error_prob\tproteinIds");
				for (int i = 0; i < NUM_ENTRIES; i++) {
					// PSMId is sequence+charge
					pw.print('A' + i);
					pw.print("+");
					pw.print(2);
					pw.print("\t");
					pw.print(0d); // q-value
					pw.print("\t");
					pw.print(0d); // PEP
					pw.print("\t");
					pw.println(i); // protein id (can't be empty)
				}
			}

			return new EncyclopediaJobData(
					Files.createTempFile(tempDir, name, ".dia").toFile(), // dia file; must exist
					new PercolatorExecutionData(
							Files.createTempFile(tempDir, name, ".features.txt").toFile(), // input tsv
							null, // fasta
							peptideOutput.toFile(), // peptide output
							Files.createTempFile(tempDir, name, ".decoys.txt").toFile(), // decoy output
							null, // protein output
							null, // protein decoy
							searchParameters
					),
					searchParameters,
					"SearchToBLIBIT",
					null,
					new EncyclopediaOneScoringFactory(searchParameters)
			);
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}

	private void assertHasPercolatorMetadata(LibraryFile file) throws IOException, SQLException {
		if (null == file) {
			throw new NullPointerException("Can't run assertions on null LibraryFile!");
		}

		final HashMap<String, String> metadata = file.getMetadata();
		assertNotNull(metadata.get(LibraryFile.PERCOLATOR_VERSION));
		assertNotNull(metadata.get("pi0"));
	}
}