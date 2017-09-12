package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TDoubleObjectHashMap;

import java.util.List;
import java.util.concurrent.BlockingQueue;

public interface PecanScoringFactory {
	public String getVersion();
	public PecanSearchParameters getParameters();
	public AbstractPecanFragmentationModel getFragmentationModel(FastaPeptideEntry peptide, AminoAcidConstants aaConstants);
	public PSMScorer getBackgroundScorer();
	public PSMPeakScorer getPecanScorer();
	public AbstractPecanScoringTask getScoringTask(PSMPeakScorer scorer, List<LibraryEntry> entries, List<Stripe> stripes, TDoubleObjectHashMap<XYPoint>[] background, PrecursorScanMap precursors, int scanAveragingMargin, BlockingQueue<PeptideScoringResult> resultsQueue);
	public PeptideScoringResultsConsumer getResultsConsumer(BlockingQueue<PeptideScoringResult> resultsQueue, StripeFileInterface diaFile);
}
