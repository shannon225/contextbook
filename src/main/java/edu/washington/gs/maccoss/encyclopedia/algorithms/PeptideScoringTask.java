package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PeptideScoringTask extends ThreadableTask<HashMap<LibraryEntry, PeptideScoringResult>> {
	/**
	 * must be immutable!
	 */
	protected final PSMScorer scorer;
	protected final ArrayList<LibraryEntry> entries;
	protected final ArrayList<FragmentScan> stripes;
	protected final PrecursorScanMap precursors;
	protected final AminoAcidConstants aaConstants;

	public PeptideScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, AminoAcidConstants aaConstants) {
		this.scorer=scorer;
		this.entries=entries;
		this.stripes=stripes;
		this.precursors=precursors;
		this.aaConstants=aaConstants;
	}

	@Override
	protected HashMap<LibraryEntry, PeptideScoringResult> process() {
		HashMap<LibraryEntry, PeptideScoringResult> map=new HashMap<LibraryEntry, PeptideScoringResult>();
		for (LibraryEntry entry : entries) {
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), aaConstants);
			TFloatFloatHashMap scoreMap=new TFloatFloatHashMap();
			
			PeptideScoringResult result=new PeptideScoringResult(entry);
			for (FragmentScan stripe : stripes) {
				float score=scorer.score(entry, stripe, predictedIsotopeDistribution, precursors);
				
				float rt=stripe.getScanStartTime();
				scoreMap.put(rt, score);
				result.addStripe(score, new float[0], stripe);
			}
			//EValueCalculator calculator=new EValueCalculator(scoreMap);
			result.setTrace(new XYTrace(scoreMap, GraphType.line, entry.getPeptideModSeq()));
			map.put(entry, result);
		}
		return map;
	}

	@Override
	public String getTaskName() {
		StringBuilder sb=new StringBuilder();
		for (LibraryEntry entry : entries) {
			if (sb.length()>0) {
				sb.append(',');
			}
			sb.append(entry.getPeptideModSeq());
		}
		return sb.toString();
	}
}
