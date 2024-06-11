package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ScoredPSM;
import edu.washington.gs.maccoss.encyclopedia.algorithms.SearchTestSupport;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;
import junit.framework.TestCase;

public class EncyclopediaTwoScoringTaskTest extends TestCase {

	private LibraryFile LIBRARY=null;
	private StripeFileInterface SINGLE_WINDOW_STRIPE_FILE=null;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		if (SINGLE_WINDOW_STRIPE_FILE==null||LIBRARY==null) {
			SINGLE_WINDOW_STRIPE_FILE = SearchTestSupport.getSingleWindowStripeFile();
			LIBRARY = SearchTestSupport.getLibrary();
		}
	}
	
	private AbstractScoringResult processPeptide(String peptideModSeq, byte precursorCharge) throws Exception {
		HashMap<String, String> paramsMap=SearchParameterParser.getDefaultParameters();
		SearchParameters parameters=SearchParameterParser.parseParameters(paramsMap);
		
		PSMScorer scorer=new EncyclopediaTwoScorer(parameters);
		ArrayList<LibraryEntry> entries=LIBRARY.getEntries(peptideModSeq, precursorCharge, true);
		
		ArrayList<FragmentScan> stripes=SINGLE_WINDOW_STRIPE_FILE.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), -Float.MAX_VALUE, Float.MAX_VALUE, true);
		Collections.sort(stripes);
		
		Range precursorIsolationRange=new Range(600.523681640625, 602.523742675781);
		float dutyCycle=1.88462436199188f;
		
		PrecursorScanMap precursors=new PrecursorScanMap(SINGLE_WINDOW_STRIPE_FILE.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		BlockingQueue<AbstractScoringResult> resultsQueue=new LinkedBlockingQueue<AbstractScoringResult>();
		
		EncyclopediaTwoScoringTask task=new EncyclopediaTwoScoringTask(scorer, entries, stripes, precursorIsolationRange, dutyCycle, precursors,
				resultsQueue, parameters);
		task.call();
		
		AbstractScoringResult result=resultsQueue.peek();
		return result;
	}
	
	public void testSinglePeptideAHWTPFEGQK() {
		try {
			AbstractScoringResult result=processPeptide("AHWTPFEGQK", (byte)2);
			assertEquals("AHWTPFEGQK", result.getEntry().getPeptideModSeq());
			assertTrue(result.getBestScore()>1f);
			ScoredPSM psm=result.getScoredMSMS();
			assertTrue(new Range(47.5f, 48.2f).contains(psm.getMSMS().getScanStartTime()/60f));
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	public void testSinglePeptideSAGFHPSGSVLAVGTVTGR() {
		try {
			AbstractScoringResult result=processPeptide("SAGFHPSGSVLAVGTVTGR", (byte)3);
			assertEquals("SAGFHPSGSVLAVGTVTGR", result.getEntry().getPeptideModSeq());
			assertTrue(result.getBestScore()>8);
			ScoredPSM psm=result.getScoredMSMS();
			assertTrue(new Range(60.6f, 61f).contains(psm.getMSMS().getScanStartTime()/60f));
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	/**
	 * Smoke Test
	 * Only a single spectrum but search the whole library
	 */
	public void testSinglePeptideSADESGQALLAAGHYASDEVREK() {
		try {
			AbstractScoringResult result=processPeptide("SADESGQALLAAGHYASDEVREK", (byte)4);
			assertEquals("SADESGQALLAAGHYASDEVREK", result.getEntry().getPeptideModSeq());
			assertTrue(result.getBestScore()>5f); // 5.66
			ScoredPSM psm=result.getScoredMSMS();
			assertTrue(new Range(58f, 58.8f).contains(psm.getMSMS().getScanStartTime()/60f)); // Not the right integration, but time is 58.3
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	/**
	 * Smoke Test
	 * Only a single spectrum but search the whole library
	 */
	public void testSinglePeptideSIYEGDESFR() {
		try {
			AbstractScoringResult result=processPeptide("SIYEGDESFR", (byte)2);
			assertEquals("SIYEGDESFR", result.getEntry().getPeptideModSeq());
			assertTrue(result.getBestScore()>4f); // 4.5910397
			ScoredPSM psm=result.getScoredMSMS();
			assertTrue(new Range(48.6f, 49.2f).contains(psm.getMSMS().getScanStartTime()/60f)); // expected time for this peptide is 49 min
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	/**
	 * Smoke Test
	 * Only a single spectrum but search the whole library
	 */
	public void testSinglePeptideRNFILDQTNVSAAAQR() {
		try {
			AbstractScoringResult result=processPeptide("RNFILDQTNVSAAAQR", (byte)3);
			assertEquals("RNFILDQTNVSAAAQR", result.getEntry().getPeptideModSeq());
			assertTrue(result.getBestScore()>8); // 8.558
			ScoredPSM psm=result.getScoredMSMS();
			assertTrue(new Range(56f, 56.5f).contains(psm.getMSMS().getScanStartTime()/60f)); // expected time for this peptide is 56.3 min
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	

	public static final String[] expected = new String[] { 
			"SSEHINEGETAMLVC[+57.021464]K",
			"NSTFSEIFKK",
			"HVSIQEAESYAESVGAK",
			"RNFILDQTNVSAAAQR",
			"VVDLMAHMASK",
			"LHLGTTQNSLTEADFR",
			"HGEVC[+57.021464]PAGWKPGSDTIKPDVQK",
			"NVIGLQMGTNR",
			"KLVIIEGDLERTEER",
			"IFC[+57.021464]C[+57.021464]HGGLSPDLQSMEQIRR",
			"LEAAIAEAEER",
			"ALQASALNAWR",
			"DLSLEEIQKK",
			"AAGPLLTDEC[+57.021464]R",
			"VAEELALEQAK",
			"IPWTAASSQLK",
			"SLEDALAEAQR",
			"EVFGSGTAC[+57.021464]QVC[+57.021464]PVHR",
			"VM[+15.994915]LGETNPADSKPGTIR",
			"GEDFPANNIVK"	
	};
	
	/**
	 * Smoke Test
	 * Only a single spectrum but search the whole library
	 */
	public void testEncyclopediaTask() {
		try {
			HashMap<String, String> paramsMap=SearchParameterParser.getDefaultParameters();
			SearchParameters parameters=SearchParameterParser.parseParameters(paramsMap);
			
			PSMScorer scorer=new EncyclopediaTwoScorer(parameters);
			ArrayList<LibraryEntry> entries=LIBRARY.getAllEntries(true, parameters.getAAConstants());
			
			ArrayList<FragmentScan> stripes=SINGLE_WINDOW_STRIPE_FILE.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);
			
			Range precursorIsolationRange=new Range(600.523681640625, 602.523742675781);
			float dutyCycle=1.88462436199188f;
			
			PrecursorScanMap precursors=new PrecursorScanMap(SINGLE_WINDOW_STRIPE_FILE.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
			BlockingQueue<AbstractScoringResult> resultsQueue=new LinkedBlockingQueue<AbstractScoringResult>();
			
			EncyclopediaTwoScoringTask task=new EncyclopediaTwoScoringTask(scorer, entries, stripes, precursorIsolationRange, dutyCycle, precursors,
					resultsQueue, parameters);

			task.call();

			TFloatArrayList scoreList=new TFloatArrayList();
			TObjectFloatHashMap<String> scoresByPeptideModSeq=new TObjectFloatHashMap<>();
			for (AbstractScoringResult result : resultsQueue) {
				LibraryEntry entry = result.getEntry();
				scoreList.add(result.getBestScore());
				scoresByPeptideModSeq.put(entry.getPeptideModSeq(), result.getBestScore());
				System.out.println(entry.getPeptideModSeq()+"\t"+result.getScoredMSMS().getMSMS().getScanStartTime()/60f+"\t"+entry.getPrecursorCharge()+"\t"+result.getBestScore());
			}
			float[] scoreArray=scoreList.toArray();
			float top5Percent=QuickMedian.select(scoreArray, 0.95f);
			//System.out.println("top 5%: "+top5Percent);
			
			for (String peptideModSeq : expected) {
				assertTrue("Problem with: "+peptideModSeq, scoresByPeptideModSeq.get(peptideModSeq)>top5Percent);
			}
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	@Override
	protected void tearDown() throws Exception {
		LIBRARY.close();
		SINGLE_WINDOW_STRIPE_FILE.close();
	}
}
