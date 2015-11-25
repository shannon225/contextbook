package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PecanLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

//@Immutable
public class PecanRawScorer implements PSMScorer {
	private final MassTolerance fragmentTolerance;
	private final AuxillaryPSMScorer auxScorer;

	public PecanRawScorer(MassTolerance fragmentTolerance, AuxillaryPSMScorer auxScorer) {
		this.fragmentTolerance = fragmentTolerance;
		this.auxScorer=auxScorer;
	}

	@Override
	public float score(LibraryEntry entry, Stripe spectrum, PrecursorScanMap precursors) {
		return General.sum(getIndividualPeakScores(entry, spectrum, true)); // dot product
	}

	@Override
	public float[] auxScore(LibraryEntry entry, Stripe spectrum, PrecursorScanMap precursors) {
		return auxScorer.score(entry, spectrum, precursors);
	}
	
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return auxScorer.getScoreNames(entry);
	}
	
	@Override
	public float[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize) {
		double[] libraryMasses=entry.getMassArray();
		float[] libraryIntensities;
		// TODO: this seems questionable that unnormalized intensities are used for individual scores while normalized intensities are used for total scores. -BCS
		if (!normalize&&entry instanceof PecanLibraryEntry) {
			libraryIntensities=((PecanLibraryEntry)entry).getUnnormalizedIntensities();
		} else {
			libraryIntensities=entry.getIntensityArray();
		}
		
		double[] spectrumMasses=spectrum.getMassArray();
		float[] spectrumIntensities=spectrum.getIntensityArray();
		
		float[] individualPeakScores=new float[libraryMasses.length];
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return individualPeakScores;
		
		for (int i=0; i<libraryMasses.length; i++) {
			int[] indicies=fragmentTolerance.getIndicies(spectrumMasses, libraryMasses[i]);
			float intensity=0.0f;
			for (int j=0; j<indicies.length; j++) {
				intensity+=spectrumIntensities[indicies[j]];
			}
			float peakScore=libraryIntensities[i]*intensity;
			individualPeakScores[i]=peakScore;
		}
		
		return individualPeakScores;
	}
}
