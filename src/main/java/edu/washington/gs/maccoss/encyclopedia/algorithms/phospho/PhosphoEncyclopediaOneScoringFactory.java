package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class PhosphoEncyclopediaOneScoringFactory implements LibraryScoringFactory {
	public static final String version="0.4.10";
	private final SearchParameters parameters;
	private final PhosphoLocalizer localizer;

	public PhosphoEncyclopediaOneScoringFactory(SearchParameters parameters, PhosphoLocalizer localizer) {
		this.parameters=parameters;
		this.localizer=localizer;
	}

	@Override
	public PSMPeakScorer getLibraryScorer(LibraryBackgroundInterface background) {
		return new EncyclopediaOneScorer(parameters, background); 
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, BlockingQueue<PeptideScoringResult> resultsQueue, File diaFile) {
		return new ScoringResultsToTSVConsumer(outputFile, diaFile, General.concatenate(EncyclopediaOneAuxillaryPSMScorer.getScoreNames(true), "localizationScore"), resultsQueue, 1);
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		return new CASiLOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, resultsQueue, parameters);
	}
	
	@Override
	public AbstractLibraryScoringTask getDDAScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		throw new EncyclopediaException("Not implemented");
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
