package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;
import junit.framework.TestCase;

public class LibraryToQuantTest extends TestCase {
	private LibraryFile RESULT_LIBRARY=null;
	private StripeFileInterface SINGLE_WINDOW_STRIPE_FILE=null;
	private SearchParameters PARAMETERS=null;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		if (SINGLE_WINDOW_STRIPE_FILE==null||RESULT_LIBRARY==null) {
			SINGLE_WINDOW_STRIPE_FILE = SearchTestSupport.getSingleWindowStripeFile();
			RESULT_LIBRARY = SearchTestSupport.getResultLibrary();
			HashMap<String, String> paramsMap=SearchParameterParser.getDefaultParameters();
			PARAMETERS=SearchParameterParser.parseParameters(paramsMap);
		}
	}
	
	public void testSaveQuantData() throws Exception {
		System.out.println(RESULT_LIBRARY.getFile().getAbsolutePath());
		System.out.println(SINGLE_WINDOW_STRIPE_FILE.getFile().getAbsolutePath());
		
		//File resultFile=new File("/Users/searleb/Downloads/precursor.elib"); 
		//File resultFile=File.createTempFile("encyclopedia_test_", ".elib");
		//resultFile.deleteOnExit();
		
		File rawDirectory=SINGLE_WINDOW_STRIPE_FILE.getFile().getParentFile();
		
		// integrate precursors
		//LibraryToQuant.saveQuantData(RESULT_LIBRARY, resultFile, rawDirectory, true, PARAMETERS, new EmptyProgressIndicator());
		Pair<TObjectFloatHashMap<String>, ArrayList<IntegratedLibraryEntry>> pair=LibraryToQuant.extractQuantData(RESULT_LIBRARY, rawDirectory, false, PARAMETERS, new EmptyProgressIndicator());
		
		// integrate fragments
		Pair<TObjectFloatHashMap<String>, ArrayList<IntegratedLibraryEntry>> pair2=LibraryToQuant.extractQuantData(RESULT_LIBRARY, rawDirectory, true, PARAMETERS, new EmptyProgressIndicator());
		
		HashMap<String, IntegratedLibraryEntry> peptideModSeqs=new HashMap<>();
		for (IntegratedLibraryEntry entry : pair.y) {
			peptideModSeqs.put(entry.getPeptideModSeq(), entry);
		}
//		for (IntegratedLibraryEntry entry : pair2.y) {
//			IntegratedLibraryEntry integratedLibraryEntry = peptideModSeqs.get(entry.getPeptideModSeq());
//			if (integratedLibraryEntry!=null) {
//				float frag = integratedLibraryEntry.getRefinementData().getQuantitativeValue();
//				float prec = entry.getRefinementData().getQuantitativeValue();
//				System.out.println(frag+","+prec);
//			}
//		}

		assertEquals(140, pair2.y.size()); // 140 of 140
		assertEquals(133, pair.y.size()); // 139 of 140
		assertTrue(peptideModSeqs.containsKey("IPWTAASSQLK")); // a correct integration
		assertFalse(peptideModSeqs.containsKey("NLENTQNQIK")); // missing because no peaks
		
		
	}
	
}
