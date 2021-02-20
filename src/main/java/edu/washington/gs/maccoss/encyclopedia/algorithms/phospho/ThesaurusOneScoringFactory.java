package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class ThesaurusOneScoringFactory implements LibraryScoringFactory {
	private final ThesaurusSearchParameters parameters;
	private final PhosphoLocalizer localizer;
	private final BlockingQueue<ModificationLocalizationData> localizationQueue;

	public ThesaurusOneScoringFactory(ThesaurusSearchParameters parameters, PhosphoLocalizer localizer, BlockingQueue<ModificationLocalizationData> localizationQueue) {
		this.parameters=parameters;
		this.localizer=localizer;
		this.localizationQueue=localizationQueue;
	}

	@Override
	public PSMPeakScorer getLibraryScorer(LibraryBackgroundInterface background) {
		return new EncyclopediaOneScorer(parameters, background); 
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, BlockingQueue<PeptideScoringResult> resultsQueue, StripeFileInterface diaFile) {
		return new ScoringResultsToTSVConsumer(outputFile, diaFile, General.concatenate(EncyclopediaOneAuxillaryPSMScorer.getScoreNames(true), "localizationScore"), resultsQueue, parameters);
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
        return new ThesaurusOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, resultsQueue, localizationQueue, parameters);
        //return new CASiLOneScoringTask(scorer, entries, stripes, dutyCycle, precursors, localizer, resultsQueue, localizationQueue, parameters);
	}
	
	@Override
	public AbstractLibraryScoringTask getDDAScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		throw new EncyclopediaException("Not implemented");
	}

	@Override
	public SearchParameters getParameters() {
		return parameters;
	}
	
	public BlockingQueue<ModificationLocalizationData> getLocalizationQueue() {
		return localizationQueue;
	}
}
