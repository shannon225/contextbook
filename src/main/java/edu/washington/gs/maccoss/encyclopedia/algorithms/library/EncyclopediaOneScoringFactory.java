package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;

public class EncyclopediaOneScoringFactory implements LibraryScoringFactory {
	public static final String version="0.4.1";
	private final SearchParameters parameters;

	public EncyclopediaOneScoringFactory(SearchParameters parameters) {
		this.parameters=parameters;
	}

	@Override
	public PSMPeakScorer getLibraryScorer(LibraryBackground background) {
		return new EncyclopediaOneScorer(parameters, background); 
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, BlockingQueue<PeptideScoringResult> resultsQueue, File diaFile) {
		return new ScoringResultsToTSVConsumer(outputFile, diaFile, EncyclopediaOneAuxillaryPSMScorer.getScoreNames(true), resultsQueue, 1);
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		return new EncyclopediaOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, resultsQueue, parameters);
	}
	
	@Override
	public AbstractLibraryScoringTask getDDAScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		return new EncyclopediaDDAScoringTask(scorer, entries, stripes, precursors, resultsQueue, parameters);
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
