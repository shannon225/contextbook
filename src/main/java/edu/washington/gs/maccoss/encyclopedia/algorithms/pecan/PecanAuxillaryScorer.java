package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import gnu.trove.list.array.TFloatArrayList;

//@Immutable
public class PecanAuxillaryScorer implements AuxillaryPSMScorer {
	private final SearchParameters parameters;
	private final float maxPPMError;

	public PecanAuxillaryScorer(SearchParameters parameters) {
		this.parameters=parameters;
		maxPPMError=(float)parameters.getPrecursorTolerance().getPpmTolerance();
	}

	@Override
	public float[] score(LibraryEntry entry, Stripe spectrum, PrecursorScanMap precursors) {
		// precursor scoring
		float[] precursorScores=getPrecursorScores(entry, spectrum.getScanStartTime(), precursors);
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
		float spectrumMagnitude=spectrum.getIntensityMagnitude();
		float peakSimilarity=spectrumMagnitude<=0?0.0f:rawScore/spectrumMagnitude; // FINAL SCORE
		
		float individualIonThreshold=rawScore/(entry.getPeptideSeq().length()+1);
		for (float peak : individualPeakScores.toArray()) {
			if (peak>individualIonThreshold) {
				numAboveThresholdMatches++;
			}
		}
		
		return new float[] {rawScore, peakSimilarity, weightedRawScore, numAboveThresholdMatches, numMatches, averageAbsPPM, averagePPM, isotopeDotProduct};
	}
	
	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		return new String[] {"rawScore", "peakSimilarity", "weightedRawScore", "numAboveThresholdMatches", "numMatches", "averageAbsPPM", "averagePPM", "isotopeDotProduct"};
	}

	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		return new float[] {0, 0, 0, 0, 0, maxPPMError, maxPPMError, 0};
	}


	public float[] getPrecursorScores(LibraryEntry entry, float spectrumRT, PrecursorScanMap precursors) {
		byte charge=entry.getPrecursorCharge();
		Peak[] precursorPacket=precursors.getIsotopePacket(entry.getPrecursorMZ(), spectrumRT, charge, parameters.getPrecursorTolerance());
		Pair<double[], float[]> pair=Peak.toArrays(precursorPacket);
		double[] masses=pair.x;
		float[] intensities=pair.y;
		
		// weighted average AbsPPM
		float averagePPM=0.0f; // FINAL SCORE
		float averageAbsPPM=0.0f; // FINAL SCORE
		float sumIntensities=0.0f;
		// start at 1 to drop "-1" isotope
		for (int i = 1; i < masses.length; i++) {
			byte isotope=(byte)(i-1);
			double predicted=entry.getPrecursorMZ()+(isotope*MassConstants.neutronMass/charge);
			
			if (intensities[i]>0) {
				double delta=masses[i]-predicted;
				float ppm=(float)((delta/entry.getPrecursorMZ())*1000000.0*intensities[i]);
				averagePPM+=ppm;
				averageAbsPPM+=Math.abs(ppm);
				sumIntensities+=intensities[i];
			}
		}
		if (sumIntensities>0) {
			averagePPM=averagePPM/sumIntensities;
			averageAbsPPM=averageAbsPPM/sumIntensities;
		} else {
			averagePPM=maxPPMError;
			averageAbsPPM=maxPPMError;
		}
		// precursor idotp
		intensities=IsotopicDistributionCalculator.normalizeToMax(intensities);
		float[] predicted=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
		float isotopeDotProduct=0.0f; // FINAL SCORE
		float euclideanDistanceIntensities=0.0f;
		float euclideanDistancePredicted=0.0f;
		for (int i = 0; i < PrecursorScanMap.isotopes.length; i++) {
			byte isotope=PrecursorScanMap.isotopes[i];
			if (isotope>=0) {
				// intensities[i] contains an extra -1 isotope
				isotopeDotProduct+=intensities[i]*predicted[isotope];
				euclideanDistanceIntensities+=intensities[i]*intensities[i];
				euclideanDistancePredicted+=predicted[isotope]*predicted[isotope];
			}
		}
		if (euclideanDistanceIntensities>0.0f&&euclideanDistancePredicted>0.0f) {
			euclideanDistanceIntensities=(float)Math.sqrt(euclideanDistanceIntensities);
			euclideanDistancePredicted=(float)Math.sqrt(euclideanDistancePredicted);
			isotopeDotProduct=isotopeDotProduct/(euclideanDistanceIntensities*euclideanDistancePredicted);
		} else {
			isotopeDotProduct=0.0f;
		}
		
		return new float[] {averageAbsPPM, isotopeDotProduct, averagePPM};
	}
}
