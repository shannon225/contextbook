package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import gnu.trove.list.array.TFloatArrayList;

//@Immutable
public class PecanAuxillaryScorer extends AuxillaryPSMScorer {

	public PecanAuxillaryScorer(SearchParameters parameters) {
		super(parameters);
	}

	@Override
	public float[] score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		// precursor scoring
		float[] precursorScores=getPrecursorScores(entry, spectrum.getScanStartTime(), predictedIsotopeDistribution, precursors);
		float averageAbsPPM=precursorScores[0]; // FINAL SCORE
		float isotopeDotProduct=precursorScores[1]; // FINAL SCORE
		float averagePPM=precursorScores[2]; // FINAL SCORE
		
		// fragment scoring
		double[] libraryMasses=entry.getMassArray();
		float[] libraryIntensities=entry.getIntensityArray();
		
		double[] spectrumMasses=spectrum.getMassArray();
		float[] spectrumIntensities=spectrum.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) {
			return getMissingDataScores(entry);
		}
		
		int numMatches=0; // FINAL SCORE
		int numAboveThresholdMatches=0; // FINAL SCORE
		float rawScore=0.0f; // FINAL SCORE
		float weightedRawScore=0.0f; // FINAL SCORE
		float sumLibraryMasses=0.0f;
		
		TFloatArrayList individualPeakScores=new TFloatArrayList();
		
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			int compare=parameters.getFragmentTolerance().compareTo(libraryMasses[libraryIndex], spectrumMasses[spectrumIndex]);
			if (compare==0) {
				numMatches++;
				float product = libraryIntensities[libraryIndex]*spectrumIntensities[spectrumIndex];
				individualPeakScores.add(product);
				
				rawScore+=product;
				weightedRawScore+=(float)product*libraryMasses[libraryIndex];
				sumLibraryMasses+=(float)libraryMasses[libraryIndex];
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
		if (sumLibraryMasses==0.0f) {
			weightedRawScore=0.0f;
		} else {
			weightedRawScore=weightedRawScore/sumLibraryMasses;
		}
		float spectrumMagnitude;
		if (spectrum instanceof Stripe) {
			spectrumMagnitude=((Stripe)spectrum).getIntensityMagnitude();
		} else {
			float magnitude=0.0f;
			for (float f : spectrum.getIntensityArray()) {
				magnitude+=f*f;
			}
			spectrumMagnitude=(float)Math.sqrt(magnitude);
		}
		
		float peakSimilarity=spectrumMagnitude<=0?0.0f:rawScore/spectrumMagnitude; // FINAL SCORE
		
		float individualIonThreshold=rawScore/(entry.getPeptideSeq().length()+1);
		for (float peak : individualPeakScores.toArray()) {
			if (peak>individualIonThreshold) {
				numAboveThresholdMatches++;
			}
		}
		
		return new float[] {rawScore, peakSimilarity, spectrumMagnitude, weightedRawScore, numAboveThresholdMatches, numMatches, averageAbsPPM, averagePPM, isotopeDotProduct};
	}
	
	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		return new String[] {"rawScore", "peakSimilarity", "weightedRawScore", "numAboveThresholdMatches", "numMatches", "averageAbsPPM", "averagePPM", "isotopeDotProduct"};
	}

	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		float maxPrePPMError=(float)parameters.getPrecursorTolerance().getToleranceThreshold();
		return new float[] {0, 0, 0, 0, 0, maxPrePPMError, maxPrePPMError, 0};
	}
}
