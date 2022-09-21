package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFile;
import edu.washington.gs.maccoss.encyclopedia.filereaders.WindowData;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;

public class EncyclopediaOneScoringTaskTest {//extends TestCase {
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
}
