package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoreCombiner;

public class CombinablePSMScorer implements PSMScorer {
	private final PSMScorer scorer;
	private final ScoreCombiner combiner;
	
	public CombinablePSMScorer(PSMScorer scorer, ScoreCombiner combiner) {
		this.scorer=scorer;
		this.combiner=combiner;
	}

	@Override
	public float[] auxScore(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return scorer.auxScore(entry, spectrum, predictedIsotopeDistribution, precursors);
	}
	
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return scorer.getAuxScoreNames(entry);
	}
	
	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize) {
		return scorer.getIndividualPeakScores(entry, spectrum, normalize);
	}
	
	@Override
	public float score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return combiner.getScore(scorer.auxScore(entry, spectrum, predictedIsotopeDistribution, precursors));
	}
}
