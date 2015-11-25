package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public interface PecanScoringFactory {
	public SearchParameters getParameters();
	public AbstractPecanFragmentationModel getFragmentationModel(String modifiedSequence, AminoAcidConstants aaConstants);
	public PSMScorer getBackgroundScorer();
	public PSMScorer getPecanScorer();
	public AbstractPecanScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, TDoubleObjectHashMap<XYPoint> background, PrecursorScanMap precursors, int scanAveragingMargin);
}
