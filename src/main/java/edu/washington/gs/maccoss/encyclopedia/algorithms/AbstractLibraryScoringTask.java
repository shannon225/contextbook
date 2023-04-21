package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ThreadableTask;

public abstract class AbstractLibraryScoringTask extends ThreadableTask<Nothing> {
	/**
	 * must be immutable!
	 */
	protected final PSMScorer scorer;
	protected final ArrayList<LibraryEntry> entries;
	protected final ArrayList<FragmentScan> stripes;
	protected final PrecursorScanMap precursors;
	protected final BlockingQueue<AbstractScoringResult> resultsQueue;
	protected final SearchParameters parameters;

	/**
	 * scorer must be a 
	 * @param scorer
	 * @param entries
	 * @param stripes
	 * @param background
	 * @param precursors
	 * @param scanAveragingMargin
	 */
	public AbstractLibraryScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue, SearchParameters parameters) {
		this.scorer=scorer;
		this.entries=entries;
		this.stripes=stripes;
		this.precursors=precursors;
		this.resultsQueue=resultsQueue;
		this.parameters=parameters;
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
