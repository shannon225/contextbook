package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.junit.AfterClass;
import org.junit.BeforeClass;

import java.io.File;

import static edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils.getResourceAsTempFile;

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
		parameters = SearchParameterParser.getDefaultParametersObject();
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

	public static SearchJobData makeAndDoJob(File dia) throws Exception {
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
	public LibraryFile[] getReferenceSearches() throws Exception {
		final File referenceSearch1File = getResourceAsTempFile(getClass(), REFERENCE_SEARCH1_RESOURCE, tempDir, "reference1_", ".elib").toFile();
		final LibraryFile referenceSearch1 = new LibraryFile();
		referenceSearch1.openFile(referenceSearch1File);

		final File referenceSearch2File = getResourceAsTempFile(getClass(), REFERENCE_SEARCH2_RESOURCE, tempDir, "reference2_", ".elib").toFile();
		final LibraryFile referenceSearch2 = new LibraryFile();
		referenceSearch2.openFile(referenceSearch2File);

		final File referenceSearch3File = getResourceAsTempFile(getClass(), REFERENCE_SEARCH3_RESOURCE, tempDir, "reference3_", ".elib").toFile();
		final LibraryFile referenceSearch3 = new LibraryFile();
		referenceSearch3.openFile(referenceSearch3File);

		return new LibraryFile[] {referenceSearch1,referenceSearch2,referenceSearch3};
	}

	@Override
	public LibraryFile getReferenceSingleQuant() throws Exception {
		final File referenceSingleQuant = getResourceAsTempFile(getClass(), REFERENCE_SINGLE_QUANT_RESOURCE, tempDir, "single_", ".elib").toFile();
		final LibraryFile libraryFile = new LibraryFile();
		libraryFile.openFile(referenceSingleQuant);
		return libraryFile;
	}

	@Override
	public LibraryFile getReferenceMulti() throws Exception {
		final File referenceMulti = getResourceAsTempFile(getClass(), REFERENCE_MULTI_RESOURCE, tempDir, "multi_", ".elib").toFile();
		final LibraryFile libraryFile = new LibraryFile();
		libraryFile.openFile(referenceMulti);
		return libraryFile;
	}

	@Override
	public LibraryFile getReferenceMultiQuant() throws Exception {
		final File referenceMultiQuant = getResourceAsTempFile(getClass(), REFERENCE_MULTI_QUANT_RESOURCE, tempDir, "multi_quant_", ".elib").toFile();
		final LibraryFile libraryFile = new LibraryFile();
		libraryFile.openFile(referenceMultiQuant);
		return libraryFile;
	}
}
