package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.pecan.PecanSearchParameters;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class VariantXCorDIAOneScoringFactory extends XCorDIAOneScoringFactory {
	public static final String version="0.1";
	private BackgroundFrequencyInterface background=null;
	private final BlockingQueue<ModificationLocalizationData> localizationQueue;

	public VariantXCorDIAOneScoringFactory(PecanSearchParameters parameters, BackgroundFrequencyInterface background, BlockingQueue<ModificationLocalizationData> localizationQueue) {
		super(parameters);
		this.background=background;
		this.localizationQueue=localizationQueue;
	}
	
	public BlockingQueue<ModificationLocalizationData> getLocalizationQueue() {
		return localizationQueue;
	}
	
	public void setBackground(BackgroundFrequencyInterface background) {
		this.background=background;
	}

	@Override
	public PSMScorer getLibraryScorer(LibraryBackgroundInterface background) {
		return new XCorDIAOneScorer(getParameters(), background); 
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, BlockingQueue<PeptideScoringResult> resultsQueue, StripeFileInterface diaFile) {
		return new ScoringResultsToTSVConsumer(outputFile, diaFile, General.concatenate(EncyclopediaOneAuxillaryPSMScorer.getScoreNames(false), "neededToLocalize", "numberOfWellShapedIons"), resultsQueue, 1);
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		if (background==null) throw new EncyclopediaException("You must initialize background before generating a scoring task!");
		return new VariantXcorDIAOneScoringTask(scorer, background, entries, stripes, precursorIsolationRange, dutyCycle, precursors, resultsQueue, localizationQueue, getParameters());
	}
	
	@Override
	public AbstractLibraryScoringTask getDDAScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue) {
		throw new EncyclopediaException("Sorry, DDA scoring for XCorDIA is not implemented!");
	}

	@Override
	public String getVersion() {
		return version;
	}
}
