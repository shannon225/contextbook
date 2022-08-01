package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class ScribeScoringFactory implements LibraryScoringFactory {
	private final SearchParameters parameters;

	public ScribeScoringFactory(SearchParameters parameters) {
		this.parameters=parameters;
	}
	
	@Override
	public String getPrimaryScoreName() {
		return EncyclopediaOneScorer.getPrimaryScoreName();
	}
	
	@Override
	public String getName() {
		return "Scribe 1.X Scoring System";
	}

	@Override
	public PSMScorer getLibraryScorer(LibraryBackgroundInterface background) {
		return new ScribeScorer(parameters); 
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, BlockingQueue<AbstractScoringResult> resultsQueue, StripeFileInterface diaFile, LibraryInterface library) {
		return new ScoringResultsToTSVConsumer(outputFile, diaFile, ScribeAuxillaryPSMScorer.getScoreNames(), resultsQueue, parameters);
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue) {
		throw new EncyclopediaException("Sorry, DIA scoring for Scribe is not implemented!");
	}
	
	@Override
	public AbstractLibraryScoringTask getDDAScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue) {
		return new ScribeScoringTask(scorer, entries, stripes, precursors, resultsQueue, parameters);
	}

	@Override
	public SearchParameters getParameters() {
		return parameters;
	}
}
