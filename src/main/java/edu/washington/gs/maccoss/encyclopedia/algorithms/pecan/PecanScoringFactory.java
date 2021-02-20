package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public interface PecanScoringFactory {
	public PecanSearchParameters getParameters();
	public AbstractPecanFragmentationModel getFragmentationModel(FastaPeptideEntry peptide, AminoAcidConstants aaConstants);
	public PSMScorer getBackgroundScorer();
	public PSMPeakScorer getPecanScorer();
	public AbstractPecanScoringTask getScoringTask(PSMPeakScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, TDoubleObjectHashMap<XYPoint>[] background, PrecursorScanMap precursors, int scanAveragingMargin, BlockingQueue<PeptideScoringResult> resultsQueue);
	public PeptideScoringResultsConsumer getResultsConsumer(BlockingQueue<PeptideScoringResult> resultsQueue, StripeFileInterface diaFile);
}
