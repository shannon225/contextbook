package edu.washington.gs.maccoss.encyclopedia;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import org.junit.After;
import org.junit.Before;

import java.io.File;

public class EncyclopediaEndToEndIT extends AbstractEndToEndIT{

	SearchParameters parameters;
	LibraryScoringFactory libraryScoringFactory;

	@Before
	public void setUp() throws Exception {
		parameters = SearchParameterParser.getDefaultParametersObject();
		libraryScoringFactory = new EncyclopediaOneScoringFactory(parameters);
		super.setUp();
	}

	@After
	public void tearDown() throws Exception {
		super.tearDown();
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
}
