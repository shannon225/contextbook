package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;

public interface PSMScorer {
	public abstract float score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors);
	public abstract float[] auxScore(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors);
	public abstract String[] getAuxScoreNames(LibraryEntry entry);
	public abstract PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize);
	public abstract PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize, double[] ions);
}