package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Swath;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PeptideScoringTask extends ThreadableTask<HashMap<LibraryEntry, XYTrace>> {
	/**
	 * must be immutable!
	 */
	private final PSMScorer scorer;
	private final ArrayList<LibraryEntry> entries;
	private final ArrayList<Swath> swaths;
	private final TDoubleObjectHashMap<XYPoint> background; // if not null, then score using zscore (otherwise use raw score)

	public PeptideScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Swath> swaths) {
		this.scorer=scorer;
		this.entries=entries;
		this.swaths=swaths;
		this.background=null;
	}
	public PeptideScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Swath> swaths, TDoubleObjectHashMap<XYPoint> background) {
		this.scorer=scorer;
		this.entries=entries;
		this.swaths=swaths;
		this.background=background;
	}

	@Override
	protected HashMap<LibraryEntry, XYTrace> process() {
		HashMap<LibraryEntry, XYTrace> map=new HashMap<LibraryEntry, XYTrace>();
		for (LibraryEntry entry : entries) {
			TFloatFloatHashMap scoreMap=new TFloatFloatHashMap();
			for (Swath swath : swaths) {
				float score=scorer.score(entry, swath);
				float rt=swath.getScanStartTime();
				if (background!=null) {
					XYPoint meanStdev=background.get((double)rt);
					if (meanStdev.y==0.0) {
						scoreMap.put(rt, 0.0f);
					} else {
						float zscore=(float)((score-meanStdev.x)/meanStdev.y);
						scoreMap.put(rt, zscore);
					}
				} else {
					scoreMap.put(rt, score);
				}
			}
			//EValueCalculator calculator=new EValueCalculator(scoreMap);
			map.put(entry, new XYTrace(scoreMap, GraphType.line, entry.getPeptideModSeq()));
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
