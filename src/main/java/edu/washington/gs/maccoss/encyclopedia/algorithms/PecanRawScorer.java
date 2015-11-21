package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import gnu.trove.list.array.TFloatArrayList;

//@Immutable
public class PecanRawScorer implements PSMScorer {
	private final MassTolerance fragmentTolerance;
	private final PecanAuxillaryScorer auxScorer;

	public PecanRawScorer(MassTolerance fragmentTolerance, PecanAuxillaryScorer auxScorer) {
		this.fragmentTolerance = fragmentTolerance;
		this.auxScorer=auxScorer;
	}
	
	public float score(LibraryEntry entry, Stripe spectrum) {
		return subScore(entry, spectrum)[0];
	}
	
	public float[] auxScore(LibraryEntry entry, Stripe spectrum) {
		return auxScorer.score(entry, spectrum);
	}
	public PecanAuxillaryScorer getAuxScorer() {
		return auxScorer;
	}
	
	public float[] subScore(LibraryEntry entry, Stripe spectrum) {
		double[] libraryMasses=entry.getMassArray();
		float[] libraryIntensities=entry.getIntensityArray();
		
		double[] spectrumMasses=spectrum.getMassArray();
		float[] spectrumIntensities=spectrum.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return new float[] {0.0f, 0.0f};

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
		
		return new float[] {rawScore, numAboveThresholdMatches};
	}
}
