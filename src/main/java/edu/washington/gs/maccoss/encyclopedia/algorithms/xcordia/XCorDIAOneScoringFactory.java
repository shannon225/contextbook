package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackground;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class XCorDIAOneScoringFactory implements LibraryScoringFactory {
	public static final String version="0.0.1";
	private final PecanSearchParameters parameters;

	public XCorDIAOneScoringFactory(PecanSearchParameters parameters) {
		this.parameters=parameters;
	}

	@Override
	public PSMScorer getLibraryScorer(LibraryBackground background) {
		return new XCorDIAOneScorer(parameters, background); 
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, BlockingQueue<PeptideScoringResult> resultsQueue, File diaFile) {
		return new ScoringResultsToTSVConsumer(outputFile, diaFile, EncyclopediaOneAuxillaryPSMScorer.getScoreNames(), resultsQueue, 1);
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		return new XCorDIAOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, resultsQueue, parameters);
	}
	
	@Override
	public AbstractLibraryScoringTask getDDAScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		throw new EncyclopediaException("Sorry, DDA scoring for XCorDIA is not implemented!");
	}

	@Override
	public SearchParameters getParameters() {
		return parameters;
	}
	
	public PecanSearchParameters getPecanParameters() {
		return parameters;
	}

	@Override
	public String getVersion() {
		return version;
	}
}
