package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoreCombiner;

public class CombinablePSMScorer implements PSMPeakScorer {
	private final PSMPeakScorer scorer;
	private final ScoreCombiner combiner;
	
	public CombinablePSMScorer(PSMPeakScorer scorer, ScoreCombiner combiner) {
		this.scorer=scorer;
		this.combiner=combiner;
	}

	@Override
	public float[] auxScore(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return scorer.auxScore(entry, spectrum, predictedIsotopeDistribution, precursors);
	}
	
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return scorer.getAuxScoreNames(entry);
	}
	
	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize) {
		return scorer.getIndividualPeakScores(entry, spectrum, normalize);
	}
	
	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize, FragmentIon[] ions) {
		return scorer.getIndividualPeakScores(entry, spectrum, normalize, ions);
	}
	
	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return combiner.getScore(scorer.auxScore(entry, spectrum, predictedIsotopeDistribution, precursors));
	}
}
