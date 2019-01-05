package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public abstract class AbstractPecanScoringTask extends AbstractLibraryScoringTask {
	/**
	 * must be immutable!
	 */
	protected final int scanAveragingWindow;
	protected final TDoubleObjectHashMap<XYPoint>[] background; // if not null, then score using zscore (otherwise use raw score)

	/**
	 * scorer must be a 
	 * @param scorer
	 * @param entries
	 * @param stripes
	 * @param background
	 * @param precursors
	 * @param scanAveragingMargin
	 */
	public AbstractPecanScoringTask(PSMPeakScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, TDoubleObjectHashMap<XYPoint>[] background, PrecursorScanMap precursors, int scanAveragingWindow, BlockingQueue<PeptideScoringResult> resultsQueue, PecanSearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.background=background;
		this.scanAveragingWindow=scanAveragingWindow;
	}
	
	PecanSearchParameters getPecanSearchParameters() {
		return (PecanSearchParameters)parameters;
	}
}
