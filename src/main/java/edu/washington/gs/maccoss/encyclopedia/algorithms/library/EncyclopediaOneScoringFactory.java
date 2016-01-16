package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;

public class EncyclopediaOneScoringFactory implements LibraryScoringFactory {
	public static final String version="0.1";
	private final SearchParameters parameters;
	private final File outputFile;

	public EncyclopediaOneScoringFactory(SearchParameters parameters, File outputFile) {
		this.parameters=parameters;
		this.outputFile=outputFile;
	}

	@Override
	public PSMScorer getLibraryScorer() {
		return new EncyclopediaOneScorer(parameters); 
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(BlockingQueue<PeptideScoringResult> resultsQueue) {
		// FIXME Auto-generated method stub
		return null;
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors,
			BlockingQueue<PeptideScoringResult> resultsQueue) {
		// FIXME Auto-generated method stub
		return null;
	}

	@Override
	public SearchParameters getParameters() {
		return parameters;
	}

	@Override
	public String getVersion() {
		return version;
	}
}
