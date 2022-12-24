package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
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
	public String getPrimaryScoreName();
	public AbstractPecanScoringTask getScoringTask(PSMPeakScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, TDoubleObjectHashMap<XYPoint>[] background, PrecursorScanMap precursors, int scanAveragingMargin, BlockingQueue<AbstractScoringResult> resultsQueue);
	public PeptideScoringResultsConsumer getResultsConsumer(BlockingQueue<AbstractScoringResult> resultsQueue, StripeFileInterface diaFile);

	public String getName();
}
