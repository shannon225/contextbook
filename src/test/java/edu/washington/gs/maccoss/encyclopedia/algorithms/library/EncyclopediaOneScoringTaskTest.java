package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
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
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;
import junit.framework.TestCase;

public class EncyclopediaOneScoringTaskTest extends TestCase {
	public static void main(String[] args) throws Exception {
		SearchParameters parameters=SearchParameterParser.getDefaultParametersObject();
		File f=new File("/Users/searleb/Documents/encyclopedia/bugs/sangtae/EXP20124_20201210XRC5_AbDIL2.wiff.dia");
		File lib=new File("/Users/searleb/Documents/encyclopedia/bugs/sangtae/Giant_PepCal_lib.dlib");
		StripeFile raw=new StripeFile(true);
		raw.openFile(f);
		
		LibraryFile library=new LibraryFile();
		library.openFile(lib);
		
		//2020dec03_cobbs_cmv_inf_gpfdia_05_3.mzML:3876.1084:NLVPMVATVQGQNLK+2
		String peptideModSeq="SGGLLWQLVR[+10.008269]";
		byte precursorCharge=2;
		double targetMz=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
		
		HashMap<Range, WindowData> ranges=raw.getRanges();
		Range range=null;
		WindowData data=null;
		for (Entry<Range, WindowData> pair : ranges.entrySet()) {
			if (pair.getKey().contains(targetMz)) {
				range=pair.getKey();
				data=pair.getValue();
				break;
			}
		}
		assert(range!=null);
		
		float dutyCycle=data.getAverageDutyCycle();
		PrecursorScanMap precursors=new PrecursorScanMap(raw.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));

		
		ArrayList<LibraryEntry> entries=library.getEntries(range, true, parameters.getAAConstants());
		LibraryBackgroundInterface background=new LibraryBackground(entries);
		
		ArrayList<LibraryEntry> targetEntries=new ArrayList<>();
		for (LibraryEntry entry : entries) {
			if (peptideModSeq.equals(entry.getPeptideModSeq())) {
				targetEntries.add(entry);
			}
		}
		assert(targetEntries.size()==1);
		System.out.println("Peptide:"+targetEntries.get(0).getPeptideModSeq());
		
		ArrayList<FragmentScan> stripes=raw.getStripes(targetMz, 0, Float.MAX_VALUE, false);

		BlockingQueue<AbstractScoringResult> resultsQueue=new LinkedBlockingQueue<AbstractScoringResult>();
		EncyclopediaOneScorer scorer=new EncyclopediaOneScorer(parameters, background);
		EncyclopediaOneScoringTask task=new EncyclopediaOneScoringTask(
				scorer, targetEntries, stripes, range, dutyCycle, precursors, resultsQueue, parameters);
		task.process();
		AbstractScoringResult poll=resultsQueue.poll();
		XYTraceInterface trace=poll.getTrace();
		
		System.out.println(trace+" --> "+(poll.getScoredMSMS().getLibraryEntry().getScanStartTime()/60)+" --> scan: "+(poll.getScoredMSMS().getMSMS().getScanStartTime()/60));
		Charter.launchChart("time", "score", true, trace);
	}

	public static final String[] expected = new String[] { "AEYTEASGPC[+57.021464]ILTPHR",
			"EVIIMATNC[+57.021464]ENC[+57.021464]GHR", "DGENVSMKDPPDLLDR", "C[+57.021464]VYTIPAHQNLVTGVK",
			"YLMDEGAHLHIYDPK", "ALQASALNAWR", "VKGDMDISLPK", "VLC[+57.021464]GGDIYVPEDPKLK",
			"LILIAC[+57.021464]GTSYHAGVATR", "EVFQIASNDHDAAINR", "EVFGSGTAC[+57.021464]QVC[+57.021464]PVHR",
			"ALQAAYGASAPSVTSAALR", "GREEWESAALQNANTK", "KYPSIIVNC[+57.021464]VEEKPK", "SAGFHPSGSVLAVGTVTGR",
			"KLVIIEGDLERTEER", "IFC[+57.021464]C[+57.021464]HGGLSPDLQSMEQIRR", "NVIGLQMGTNR", "NRPSSGSLIQVVTTEGR",
			"VVDLMAHMASK", "VM[+15.994915]LGETNPADSKPGTIR", "LHLGTTQNSLTEADFR", "HVSIQEAESYAESVGAK",
			"SSEHINEGETAMLVC[+57.021464]K", "FAAATGATPIAGR", "AAGFKDPLLASGTDGVGTK", "RNFILDQTNVSAAAQR",
			"HGEVC[+57.021464]PAGWKPGSETIIPDPAGK", "LNEAKEEFTSGGPLGQK", "HGEVC[+57.021464]PAGWKPGSDTIKPDVQK" };
	
	/**
	 * Smoke Test
	 * Only a single spectrum but search the whole library
	 */
	public void testEncyclopediaTask() {
		try {
			HashMap<String, String> paramsMap=SearchParameterParser.getDefaultParameters();
			SearchParameters parameters=SearchParameterParser.parseParameters(paramsMap);
			
			PSMScorer scorer=new EncyclopediaTwoScorer(parameters);
			ArrayList<LibraryEntry> entries=SearchTestSupport.getLibrary().getAllEntries(true, parameters.getAAConstants());
	
			StripeFileInterface singleWindowStripeFile = SearchTestSupport.getSingleWindowStripeFile();
			
			ArrayList<FragmentScan> stripes=singleWindowStripeFile.getStripes(new Range(-Float.MAX_VALUE, Float.MAX_VALUE), -Float.MAX_VALUE, Float.MAX_VALUE, true);
			Collections.sort(stripes);
			
			Range precursorIsolationRange=new Range(600.523681640625, 602.523742675781);
			float dutyCycle=1.88462436199188f;
			
			PrecursorScanMap precursors=new PrecursorScanMap(singleWindowStripeFile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
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
			}
			float[] scoreArray=scoreList.toArray();
			float top5Percent=QuickMedian.select(scoreArray, 0.95f);
			
			for (String peptideModSeq : expected) {
				assertTrue(scoresByPeptideModSeq.get(peptideModSeq)>top5Percent);
			}
		} catch (Exception e) {
			fail(e.getMessage());
		}
	}
}
