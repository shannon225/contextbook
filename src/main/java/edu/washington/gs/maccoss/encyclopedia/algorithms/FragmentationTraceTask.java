package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanRawScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class FragmentationTraceTask extends PeptideScoringTask {

	public FragmentationTraceTask(PecanRawScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors) {
		super(scorer, entries, stripes, precursors);
	}
	private PecanRawScorer getScorer() {
		return (PecanRawScorer)super.scorer;
	}

	@Override
	protected HashMap<LibraryEntry, PeptideScoringResult> process() {
		HashMap<LibraryEntry, PeptideScoringResult> map=new HashMap<LibraryEntry, PeptideScoringResult>();
		
		for (LibraryEntry entry : super.entries) {
			String[] scoreNames=getScorer().getAuxScoreNames(entry);
			
			@SuppressWarnings("unchecked")
			ArrayList<XYPoint>[] dataPoints=new ArrayList[scoreNames.length];
			for (int i=0; i<dataPoints.length; i++) {
				dataPoints[i]=new ArrayList<XYPoint>();
			}
			for (Stripe stripe : super.stripes) {
				float[] fragmentScores=getScorer().auxScore(entry, stripe, precursors);
				float rt=stripe.getScanStartTime();
				for (int i=0; i<fragmentScores.length; i++) {
					dataPoints[i].add(new XYPoint(rt/60.0f, fragmentScores[i]));
				}
			}
			
			FragmentationScoringResult result=new FragmentationScoringResult(entry);
			for (int i=0; i<dataPoints.length; i++) {
				XYTrace trace=new XYTrace(dataPoints[i], GraphType.line, scoreNames[i]);
				result.addFragmentationTrace(trace);
			}

			map.put(entry, result);
		}
		return map;
	}
}
