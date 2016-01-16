package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;

public class EncyclopediaOneScoringTask extends AbstractLibraryScoringTask {

	public EncyclopediaOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue,
			SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
	}

	@Override
	protected Nothing process() {
		// TODO Auto-generated method stub
		return null;
	}

	
}
