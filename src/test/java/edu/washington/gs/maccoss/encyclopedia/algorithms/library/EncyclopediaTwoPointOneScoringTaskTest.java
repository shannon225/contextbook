package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.DataFormatException;

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
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import junit.framework.TestCase;

public class EncyclopediaTwoPointOneScoringTaskTest extends TestCase {

	private LibraryFile LIBRARY=null;
	private StripeFileInterface SINGLE_WINDOW_STRIPE_FILE=null;
	
	public static void main(String[] args) throws Exception {
		AbstractScoringResult result = testSearch("VVDESDETENQEEKAK", (byte)3);
		
		System.out.println(result.getBestScore());
	}

	private static AbstractScoringResult testSearch(String peptideModSeq, byte precursorCharge) throws IOException, SQLException, DataFormatException {
		File diaFile=new File("/Users/searleb/Documents/encyclopedia/astral/20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1.dia");
		File libFile=new File("/Users/searleb/Documents/encyclopedia/astral/HeLa.dlib");
		LibraryFile library=new LibraryFile();
		library.openFile(libFile);
		
		StripeFileInterface dia=new StripeFile();
		dia.openFile(diaFile);
		
		HashMap<String, String> paramsMap=SearchParameterParser.getDefaultParameters();
		SearchParameters parameters=SearchParameterParser.parseParameters(paramsMap);
		
		double pepMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		
		PSMScorer scorer=new EncyclopediaTwoScorer(parameters);
		ArrayList<LibraryEntry> targets=library.getEntries(peptideModSeq, precursorCharge, true);
		
		ArrayList<LibraryEntry> tasks=new ArrayList<LibraryEntry>();
		for (LibraryEntry entry : targets) {
			tasks.add(entry);
			tasks.add(entry.getDecoy(parameters));
		}
		
		ArrayList<FragmentScan> stripes=dia.getStripes(pepMz, -Float.MAX_VALUE, Float.MAX_VALUE, true);
		Collections.sort(stripes);
		
		Range precursorIsolationRange=new Range(600.523681640625, 602.523742675781);
		float dutyCycle=1.88462436199188f;
		
		PrecursorScanMap precursors=new PrecursorScanMap(dia.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		BlockingQueue<AbstractScoringResult> resultsQueue=new LinkedBlockingQueue<AbstractScoringResult>();
		
		EncyclopediaTwoPointOneScoringTask task=new EncyclopediaTwoPointOneScoringTask(scorer, tasks, stripes, precursorIsolationRange, dutyCycle, precursors,
				resultsQueue, parameters);
		task.call();
		
		AbstractScoringResult result=resultsQueue.peek();
		return result;
	}

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
		
		EncyclopediaTwoPointOneScoringTask task=new EncyclopediaTwoPointOneScoringTask(scorer, entries, stripes, precursorIsolationRange, dutyCycle, precursors,
				resultsQueue, parameters);
		task.call();
		
		AbstractScoringResult result=resultsQueue.peek();
		return result;
	}
	
	public void testSinglePeptideAHWTPFEGQK() {
		try {
			AbstractScoringResult result=processPeptide("AHWTPFEGQK", (byte)2);
			assertTrue(result.getBestScore()<1f);
			assertTrue(result.hasScoredResults());
			ScoredPSM psm=result.getScoredMSMS();
			assertNotNull(psm);
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	public void testSinglePeptideNSTFSEIFKK() {
		try {
			AbstractScoringResult result=processPeptide("NSTFSEIFKK", (byte)2);
			assertEquals("NSTFSEIFKK", result.getEntry().getPeptideModSeq());
			assertTrue(result.getBestScore()>10f);
			ScoredPSM psm=result.getScoredMSMS();
			assertTrue(new Range(56f, 57f).contains(psm.getMSMS().getScanStartTime()/60f));
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	public void testSinglePeptideSSEHINEGETAMLVCK() {
		try {
			AbstractScoringResult result=processPeptide("SSEHINEGETAMLVC[+57.021464]K", (byte)3);
			assertEquals("SSEHINEGETAMLVC[+57.021464]K", result.getEntry().getPeptideModSeq());
			assertTrue(result.getBestScore()>10f);
			ScoredPSM psm=result.getScoredMSMS();
			assertTrue(new Range(49f, 50f).contains(psm.getMSMS().getScanStartTime()/60f));
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	public void testSinglePeptideSAGFHPSGSVLAVGTVTGR() {
		try {
			AbstractScoringResult result=processPeptide("SAGFHPSGSVLAVGTVTGR", (byte)3);
			assertEquals("SAGFHPSGSVLAVGTVTGR", result.getEntry().getPeptideModSeq());
			assertTrue(result.getBestScore()>10f);
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
			assertTrue(result.getBestScore()>5f); // 6.475733
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
			assertTrue(result.getBestScore()>10f); // 14.534062
			ScoredPSM psm=result.getScoredMSMS();
			assertTrue(new Range(56f, 56.5f).contains(psm.getMSMS().getScanStartTime()/60f)); // expected time for this peptide is 56.3 min
			
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	public static final String expected = "SSEHINEGETAMLVC[+57.021464]K";
	
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
			
			EncyclopediaTwoPointOneScoringTask task=new EncyclopediaTwoPointOneScoringTask(scorer, entries, stripes, precursorIsolationRange, dutyCycle, precursors,
					resultsQueue, parameters);

			task.call();
			
			// only the top match is saved (TDC)
			assertTrue(resultsQueue.size()==1);
			assertEquals(resultsQueue.take().getEntry().getPeptideModSeq(), "SSEHINEGETAMLVC[+57.021464]K");

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
