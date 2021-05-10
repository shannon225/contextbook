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

	static SearchParameters parameters;
	static LibraryScoringFactory libraryScoringFactory;

	@BeforeClass
	public static void setUpParameters() throws Exception {
		parameters = SearchParameterParser.getDefaultParametersObject();
		libraryScoringFactory = new EncyclopediaOneScoringFactory(parameters);
		buildReports();
	}

	@AfterClass
	public static void tearDownParameters() throws Exception {
		tearDownReports();
		parameters = null;
		libraryScoringFactory = null;
	}

	@Override
	public SearchJobData makeAndDoJob(File dia) throws Exception {
		EncyclopediaJobData jobData = new EncyclopediaJobData(dia,fastaFile,libraryInterface,libraryScoringFactory);
		Encyclopedia.runSearch(new EmptyProgressIndicator(),jobData);
		return jobData;
	}

	@Override
	public int getPeptideFloor() {
		return 400;
	}

	@Override
	public int getProteinFloor() {
		return 300;
	}

	@Override
	public LibraryFile[] getReferenceSearches() throws Exception {
		File referenceSearch1File = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/encyc_dia_1.elib", tempDir, "reference1_", ".elib").toFile();
		LibraryFile referenceSearch1 = new LibraryFile() {{openFile(referenceSearch1File);}};

		File referenceSearch2File = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/encyc_dia_2.elib", tempDir, "reference2_", ".elib").toFile();
		LibraryFile referenceSearch2 = new LibraryFile() {{openFile(referenceSearch2File);}};

		File referenceSearch3File = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/encyc_dia_3.elib", tempDir, "reference3_", ".elib").toFile();
		LibraryFile referenceSearch3 = new LibraryFile() {{openFile(referenceSearch3File);}};

		return new LibraryFile[] {referenceSearch1,referenceSearch2,referenceSearch3};
	}

	@Override
	public LibraryFile getReferenceSingleQuant() throws Exception {
		File referenceSingleQuant = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/encyc_single.elib", tempDir, "single_", ".elib").toFile();
		return new LibraryFile() {{openFile(referenceSingleQuant);}};
	}

	@Override
	public LibraryFile getReferenceMulti() throws Exception {
		File referenceMulti = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/encyc_multi.elib", tempDir, "multi_", ".elib").toFile();
		return new LibraryFile() {{openFile(referenceMulti);}};
	}

	@Override
	public LibraryFile getReferenceMultiQuant() throws Exception {
		File referenceMultiQuant = getResourceAsTempFile(getClass(), "/edu/washington/gs/maccoss/encyclopedia/current_reports/encyc_multi_quant.elib", tempDir, "multi_quant_", ".elib").toFile();
		return new LibraryFile() {{openFile(referenceMultiQuant);}};
	}
}
