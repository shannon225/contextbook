package edu.washington.gs.maccoss.encyclopedia;

import com.google.common.collect.ImmutableList;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.tests.EncyclopediaTestUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.junit.AfterClass;
import org.junit.BeforeClass;

import java.io.File;
import java.nio.file.Path;
import java.util.List;

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

		final LibraryFile reference = new LibraryFile();
		try {
			final Path refElib = EncyclopediaTestUtils.getResourceAsTempFile(EncyclopediaEndToEndIT.class, REFERENCE_SEARCH1_RESOURCE, tempDir, "reference_", ".elib");
			reference.openFile(refElib.toFile());
			parameters = SearchParameterParser.parseParameters(reference.getMetadata());
		} finally {
			reference.close();
		}

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
}
