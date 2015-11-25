package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public class PecanOneScoringFactory implements PecanScoringFactory {
	private final SearchParameters parameters;

	public PecanOneScoringFactory(SearchParameters parameters) {
		this.parameters=parameters;
	}
	
	@Override
	public SearchParameters getParameters() {
		return parameters;
	}
	
	@Override
	public AbstractPecanFragmentationModel getFragmentationModel(String modifiedSequence, AminoAcidConstants aaConstants) {
		return new PecanOneFragmentationModel(modifiedSequence, aaConstants);
	}

	@Override
	public PSMScorer getBackgroundScorer() {
		return new DotProduct(parameters.getFragmentTolerance());
	}

	@Override
	public PSMScorer getPecanScorer() {
		return new PecanRawScorer(parameters.getFragmentTolerance(), new PecanAuxillaryScorer(parameters));
	}

	@Override
	public AbstractPecanScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, TDoubleObjectHashMap<XYPoint> background, PrecursorScanMap precursors,
			int scanAveragingMargin) {
		return new PecanOneScoringTask(scorer, entries, stripes, background, precursors, scanAveragingMargin);
	}

}
