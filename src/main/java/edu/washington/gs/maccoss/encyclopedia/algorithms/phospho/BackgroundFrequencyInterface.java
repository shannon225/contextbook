package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

public interface BackgroundFrequencyInterface {

	/**
	 * for display only
	 * @param precursorMz
	 * @return
	 */
	Pair<double[], float[]> getRoundedMassCounters(double precursorMz, MassTolerance tolerance);

	float[] getFrequencies(double[] ions, double precursorMz, MassTolerance tolerance);

	public LibraryBackgroundInterface getLibraryBackground(double precursorMz, final MassTolerance tolerance);
}