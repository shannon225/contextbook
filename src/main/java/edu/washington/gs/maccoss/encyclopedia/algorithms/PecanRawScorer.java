package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import gnu.trove.list.array.TFloatArrayList;

//@Immutable
public class PecanRawScorer implements PSMScorer {
	private final MassTolerance fragmentTolerance;

	public PecanRawScorer(MassTolerance fragmentTolerance) {
		this.fragmentTolerance = fragmentTolerance;
	}
	public float score(LibraryEntry entry, Stripe spectrum) {
		double[] libraryMasses=entry.getMassArray();
		float[] libraryIntensities=entry.getIntensityArray();
		
		double[] spectrumMasses=spectrum.getMassArray();
		float[] spectrumIntensities=spectrum.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return 0.0f;

		int numMatches=0; // FINAL SCORE
		int numAboveThresholdMatches=0; // FINAL SCORE
		float rawScore=0.0f; // FINAL SCORE
		
		TFloatArrayList individualPeakScores=new TFloatArrayList();
		
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			int compare=fragmentTolerance.compareTo(libraryMasses[libraryIndex], spectrumMasses[spectrumIndex]);
			if (compare==0) {
				float peakScore=libraryIntensities[libraryIndex]*spectrumIntensities[spectrumIndex];
				individualPeakScores.add(peakScore);
				rawScore+=peakScore;
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
		
		int peptideLength=entry.getPeptideSeq().length();
		float individualIonThreshold=rawScore/(peptideLength+1);
		for (float peak : individualPeakScores.toArray()) {
			if (peak>individualIonThreshold) {
				numAboveThresholdMatches++;
			}
		}
		
		if (numAboveThresholdMatches<0.5*peptideLength) {
			rawScore=0.0f;
		}
		
		return rawScore;
	}
}
