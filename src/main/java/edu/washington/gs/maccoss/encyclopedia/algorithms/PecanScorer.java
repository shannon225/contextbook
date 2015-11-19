package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import gnu.trove.list.array.TFloatArrayList;

//@Immutable
public class PecanScorer implements PSMScorer {
	private final MassTolerance fragmentTolerance;
	private final MassTolerance precursorTolerance;
	private final PrecursorScanMap precursors;


	public PecanScorer(MassTolerance fragmentTolerance, MassTolerance precursorTolerance, PrecursorScanMap precursors) {
		this.fragmentTolerance=fragmentTolerance;
		this.precursorTolerance=precursorTolerance;
		this.precursors=precursors;
	}


	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer#score(edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry, edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe)
	 */
	public float score(LibraryEntry entry, Stripe spectrum) {
		// precursor scoring
		float[] precursorScores=getPrecursorScores(entry, spectrum.getScanStartTime());
		float averageAbsPPM=precursorScores[0]; // FINAL SCORE
		float isotopeDotProduct=precursorScores[1]; // FINAL SCORE
		
		// fragment scoring
		double[] libraryMasses=entry.getMassArray();
		float[] libraryIntensities=entry.getIntensityArray();
		
		double[] spectrumMasses=spectrum.getMassArray();
		float[] spectrumIntensities=spectrum.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return 0.0f;
		
		int numMatches=0; // FINAL SCORE
		int numAboveThresholdMatches=0; // FINAL SCORE
		float rawScore=0.0f; // FINAL SCORE
		float weightedRawScore=0.0f; // FINAL SCORE
		float sumLibraryMasses=0.0f;
		
		TFloatArrayList individualPeakScores=new TFloatArrayList();
		
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			int compare=fragmentTolerance.compareTo(libraryMasses[libraryIndex], spectrumMasses[spectrumIndex]);
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
		weightedRawScore=weightedRawScore/sumLibraryMasses;
		float peakSimilarity=rawScore/spectrum.getIntensityMagnitude(); // FINAL SCORE
		
		float individualIonThreshold=rawScore/(entry.getPeptideSeq().length()+1);
		for (float peak : individualPeakScores.toArray()) {
			if (peak>individualIonThreshold) {
				numAboveThresholdMatches++;
			}
		}
		
		return rawScore;
	}


	public float[] getPrecursorScores(LibraryEntry entry, float spectrumRT) {
		Peak[] precursorPacket=precursors.getIsotopePacket(entry.getPrecursorMZ(), spectrumRT, entry.getPrecursorCharge(), precursorTolerance);
		Pair<double[], float[]> pair=Peak.toArrays(precursorPacket);
		double[] masses=pair.x;
		float[] intensities=pair.y;
		
		// weighted average AbsPPM
		float averageAbsPPM=0.0f; // FINAL SCORE
		float sumIntensities=0.0f;
		for (int i = 0; i < masses.length; i++) {
			averageAbsPPM+=(float)Math.abs((masses[i]-entry.getPrecursorMZ())/1000000.0f)*intensities[i];
			sumIntensities+=intensities[i];
		}
		averageAbsPPM=averageAbsPPM/sumIntensities;
		
		// precursor idotp
		intensities=IsotopicDistributionCalculator.normalizeToMax(intensities);
		float[] predicted=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq());
		float isotopeDotProduct=0.0f; // FINAL SCORE
		for (int i = 0; i < PrecursorScanMap.isotopes.length; i++) {
			byte isotope=PrecursorScanMap.isotopes[i];
			if (isotope>=0) {
				// intensities[i] contains an extra -1 isotope
				isotopeDotProduct+=intensities[i]*predicted[isotope];
			}
		}
		
		return new float[] {averageAbsPPM, isotopeDotProduct};
	}
}
