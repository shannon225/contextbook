package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.XCorDIAOneScoringFactory;
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

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.Assert.*;

public class SearchToBLIBIT {
	private final ProgressIndicator progress = new EmptyProgressIndicator();

	private SearchParameters searchParameters;
	private Path tempDir;

	@Before
	public void setUp() throws Exception {
		searchParameters = SearchParameterParser.getDefaultParametersObject();

		tempDir = Files.createTempDirectory("SearchToBLIBIT_");
		FileUtils.forceDeleteOnExit(tempDir.toFile());
	}

	@After
	public void tearDown() throws Exception {
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

		fail("TODO");
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

		fail("TODO");
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
		//TODO: move to resources
		final Path dia = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/LongSwath_UPS1_40fm_Ecoli_1ug-rep1.dia");
		final Path featuresTxt = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep1.dia.features.txt");
		final Path fasta = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/ups-protein-standards.fasta");
		final Path peptideOutput = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep1.dia.xcordia.txt");
		final Path decoyOutput = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep1.dia.xcordia.decoy.txt");

		return makeXCorDIAJobData(dia, featuresTxt, fasta, peptideOutput, decoyOutput);
	}

	private SearchJobData getSearchJobDataB() throws IOException, SQLException {
		//TODO: move to resources
		final Path dia = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/LongSwath_UPS1_40fm_Ecoli_1ug-rep2.dia");
		final Path featuresTxt = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep2.dia.features.txt");
		final Path fasta = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/ups-protein-standards.fasta");
		final Path peptideOutput = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep2.dia.xcordia.txt");
		final Path decoyOutput = Paths.get("/media/data/sethjust/proteomesoft/xcordia-testing/process_run_20201207-0/LongSwath_UPS1_40fm_Ecoli_1ug-rep2.dia.xcordia.decoy.txt");

		return makeXCorDIAJobData(dia, featuresTxt, fasta, peptideOutput, decoyOutput);
	}

	private XCorDIAJobData makeXCorDIAJobData(Path dia, Path featuresTxt, Path fasta, Path peptideOutput, Path decoyOutput) throws IOException, SQLException {
		Assume.assumeTrue(Files.exists(dia));

		final StripeFile diaReader = new StripeFile(true) ;
		diaReader.openFile(dia.toFile());

		return new XCorDIAJobData(
				Optional.empty(),
				Optional.empty(),
				dia.toFile(), // dia file; must exist
				diaReader,
				fasta.toFile(),
				new PercolatorExecutionData(
						featuresTxt.toFile(), // input tsv
						fasta.toFile(), // fasta
						peptideOutput.toFile(), // peptide output
						decoyOutput.toFile(), // decoy output
						null, // protein output
						null, // protein decoy
						searchParameters
				),
				new XCorDIAOneScoringFactory(new PecanSearchParameters(
						searchParameters.getAAConstants(),
						searchParameters.getFragType(),
						searchParameters.getFragmentTolerance(),
						searchParameters.getPrecursorTolerance(),
						searchParameters.getEnzyme(),
						1,
						searchParameters.isQuantifySameFragmentsAcrossSamples(),
						false,
						false
				))
		);
	}
}