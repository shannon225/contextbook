package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

//@Immutable
public class DotProduct implements PSMScorer {
	private final MassTolerance tolerance;

	public DotProduct(MassTolerance tolerance) {
		this.tolerance = tolerance;
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer#score(edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry, edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe)
	 */
	public float score(LibraryEntry entry, Stripe spectrum, PrecursorScanMap precursors) {
		float[] peakscores=getIndividualPeakScores(entry, spectrum, false);
		return General.sum(peakscores);
	}
	
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return new String[0];
	}
	@Override
	public float[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize) {
		double[] libraryMasses=entry.getMassArray();
		float[] libraryIntensities=entry.getIntensityArray();
		
		double[] spectrumMasses=spectrum.getMassArray();
		float[] spectrumIntensities=spectrum.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return new float[0];
		
		float[] peakscores=new float[libraryIntensities.length];
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			int compare=tolerance.compareTo(libraryMasses[libraryIndex], spectrumMasses[spectrumIndex]);
			if (compare==0) {
				peakscores[libraryIndex]=libraryIntensities[libraryIndex]*spectrumIntensities[spectrumIndex];
				libraryIndex++;
				spectrumIndex++;
			} else if (compare>0) {
				spectrumIndex++;
			} else {
				libraryIndex++;
			}
			if (libraryIndex>=libraryMasses.length) break;
			if (spectrumIndex>=spectrumMasses.length) break;
		}
		
		return peakscores;
	}
	
	@Override
	public float[] auxScore(LibraryEntry entry, Stripe spectrum, PrecursorScanMap precursors) {
		return new float[0];
	}
}
