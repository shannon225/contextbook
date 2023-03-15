package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ScoredPSM;
import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import junit.framework.TestCase;

public class ScribeScoringTaskTest extends TestCase {

	/**
	 * Smoke Test
	 * Only a single spectrum but search the whole library
	 */
	public void testScribeTask() {
		try {
			HashMap<String, String> paramsMap=SearchParameterParser.getDefaultParameters();
			paramsMap.put("-ptolunits", "AMU");
			SearchParameters parameters=SearchParameterParser.parseParameters(paramsMap);
			
			PSMScorer scorer=new ScribeScorer(parameters);
			ArrayList<LibraryEntry> entries=SearchTestSupport.getLibrary().getAllEntries(true, parameters.getAAConstants());
	
			StripeFileInterface singleWindowStripeFile = SearchTestSupport.getSingleWindowStripeFile();
			//ArrayList<FragmentScan> stripes=singleWindowStripeFile.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			ArrayList<FragmentScan> stripes=singleWindowStripeFile.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), 3180.3115f-1f, 3180.3115f+1f, true);
			Collections.sort(stripes);
			
			PrecursorScanMap precursors=new PrecursorScanMap(singleWindowStripeFile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
			BlockingQueue<AbstractScoringResult> resultsQueue=new LinkedBlockingQueue<AbstractScoringResult>();
			
			ScribeScoringTask task=new ScribeScoringTask(scorer, entries, stripes, precursors, resultsQueue, parameters);

			task.call();

			assertEquals("IQAVIDAGVC[+57.021464]R", resultsQueue.peek().getEntry().getPeptideModSeq());
		} catch (Exception e) {
			fail(e.getMessage());
		}
	}

}
