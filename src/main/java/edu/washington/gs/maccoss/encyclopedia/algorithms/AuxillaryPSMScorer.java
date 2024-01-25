package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public abstract class AuxillaryPSMScorer {
	public static final int MISSING_INDEX=-1;
	protected final SearchParameters parameters;
	private final float maxPPMError;

	public AuxillaryPSMScorer(SearchParameters parameters) {
		this.parameters=parameters;
		maxPPMError=(float)parameters.getPrecursorTolerance().getToleranceThreshold();
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer#score(edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry, edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe)
	 */
	public abstract float[] score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors);
	public abstract float[] getMissingDataScores(LibraryEntry entry);
	public abstract String[] getScoreNames(LibraryEntry entry);

	public float[] getPrecursorScores(LibraryEntry entry, float spectrumRT, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		byte charge=entry.getPrecursorCharge();
		Peak[] precursorPacket=precursors.getIsotopePacket(entry.getPrecursorMZ(), spectrumRT, charge, parameters.getPrecursorTolerance());
		if (precursorPacket.length==0) {
			return new float[] {maxPPMError, 0.0f, maxPPMError};
		}
		
		Triplet<double[], float[], Optional<float[]>> pair=Peak.toArrays(precursorPacket);
		double[] masses=pair.x;
		float[] intensities=pair.y;
		
		// weighted average AbsPPM
		float averagePPM=0.0f; // FINAL SCORE
		float averageAbsPPM=0.0f; // FINAL SCORE
		int peaksUsed=0;
		// start at 1 to drop "-1" isotope
		for (int i = 1; i < masses.length; i++) {
			byte isotope=(byte)(i-1);
			double predicted=entry.getPrecursorMZ()+(isotope*MassConstants.neutronMass/charge);
			
			if (intensities[i]>0) {
				// FIXME THIS SHOULD BE WEIGHTED BY INTENSITY!
				float ppm=(float)parameters.getPrecursorTolerance().getDeltaScore(predicted, masses[i]);
				averagePPM+=ppm;
				averageAbsPPM+=Math.abs(ppm);
				peaksUsed++;
			}
		}
		if (peaksUsed>0) {
			averagePPM=averagePPM/peaksUsed;
			averageAbsPPM=averageAbsPPM/peaksUsed;
		} else {
			averagePPM=maxPPMError;
			averageAbsPPM=maxPPMError;
		}
		
		int numberMatch=0;
		float percentBlankOverMono=intensities[1]<=0?1.0f:Math.min(1.0f, intensities[0]/intensities[1]);
		// precursor idotp
		intensities=IsotopicDistributionCalculator.normalizeToMax(intensities);
		float isotopeDotProduct=0.0f; // FINAL SCORE
		float euclideanDistanceIntensities=0.0f;
		float euclideanDistancePredicted=0.0f;
		for (int i = 0; i < PrecursorScanMap.isotopes.length; i++) {
			byte isotope=PrecursorScanMap.isotopes[i];
			euclideanDistanceIntensities+=intensities[i]*intensities[i]; // add -1 into iDotP
			
			if (isotope>=0) {
				if (intensities[i]>0) {
					numberMatch++;
				}
				// intensities[i] contains an extra -1 isotope
				isotopeDotProduct+=intensities[i]*predictedIsotopeDistribution[isotope];
				euclideanDistancePredicted+=predictedIsotopeDistribution[isotope]*predictedIsotopeDistribution[isotope];
			}
		}
		if (euclideanDistanceIntensities>0.0f&&euclideanDistancePredicted>0.0f) {
			euclideanDistanceIntensities=(float)Math.sqrt(euclideanDistanceIntensities);
			euclideanDistancePredicted=(float)Math.sqrt(euclideanDistancePredicted);
			isotopeDotProduct=isotopeDotProduct/(euclideanDistanceIntensities*euclideanDistancePredicted);
		} else {
			isotopeDotProduct=0.0f;
		}

		return new float[] {averageAbsPPM, isotopeDotProduct, averagePPM, percentBlankOverMono, numberMatch};
	}


	public abstract int getParentDeltaMassIndex();
	public abstract int getFragmentDeltaMassIndex();
}