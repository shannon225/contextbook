package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;

public abstract class AuxillaryPSMScorer {
	protected final SearchParameters parameters;
	private final float maxPPMError;

	public AuxillaryPSMScorer(SearchParameters parameters) {
		this.parameters=parameters;
		maxPPMError=(float)parameters.getPrecursorTolerance().getPpmTolerance();
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer#score(edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry, edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe)
	 */
	public abstract float[] score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors);
	public abstract float[] getMissingDataScores(LibraryEntry entry);
	public abstract String[] getScoreNames(LibraryEntry entry);

	public float[] getPrecursorScores(LibraryEntry entry, float spectrumRT, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		byte charge=entry.getPrecursorCharge();
		Peak[] precursorPacket=precursors.getIsotopePacket(entry.getPrecursorMZ(), spectrumRT, charge, parameters.getPrecursorTolerance());
		if (precursorPacket.length==0) {
			return new float[] {maxPPMError, 0.0f, maxPPMError};
		}
		
		Pair<double[], float[]> pair=Peak.toArrays(precursorPacket);
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
				double delta=predicted-masses[i];
				float ppm=(float)((delta/entry.getPrecursorMZ())*1000000.0);
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
		
		// precursor idotp
		intensities=IsotopicDistributionCalculator.normalizeToMax(intensities);
		float isotopeDotProduct=0.0f; // FINAL SCORE
		float euclideanDistanceIntensities=0.0f;
		float euclideanDistancePredicted=0.0f;
		for (int i = 0; i < PrecursorScanMap.isotopes.length; i++) {
			byte isotope=PrecursorScanMap.isotopes[i];
			if (isotope>=0) {
				// intensities[i] contains an extra -1 isotope
				isotopeDotProduct+=intensities[i]*predictedIsotopeDistribution[isotope];
				euclideanDistanceIntensities+=intensities[i]*intensities[i];
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

		return new float[] {averageAbsPPM, isotopeDotProduct, averagePPM};
	}
}