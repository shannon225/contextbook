package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IonType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;

//@Immutable
public class DotProduct implements PSMScorer {
	private final MassTolerance tolerance;

	public DotProduct(MassTolerance tolerance) {
		this.tolerance = tolerance;
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer#score(edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry, edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe)
	 */
	public float score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return PeakScores.sumScores(getIndividualPeakScores(entry, spectrum, false));
	}
	
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return new String[0];
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize) {
		return getIndividualPeakScores(entry, spectrum, normalize, null);
	}
	
	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize, FragmentIon[] ions) {
		if (ions!=null) {
			throw new EncyclopediaException("DotProduct doesn't currently handle ion selection. Please report this bug!");
		}
		
		double[] libraryMasses=entry.getMassArray();
		float[] libraryIntensities=entry.getIntensityArray();
		
		double[] spectrumMasses=spectrum.getMassArray();
		float[] spectrumIntensities=spectrum.getIntensityArray();
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return new PeakScores[0];
		
		PeakScores[] peakscores=new PeakScores[libraryIntensities.length];
		int libraryIndex=0;
		int spectrumIndex=0;
		while (true) {
			double targetMass=libraryMasses[libraryIndex];
			int compare=tolerance.compareTo(targetMass, spectrumMasses[spectrumIndex]);
			if (compare==0) {
				float score=libraryIntensities[libraryIndex]*spectrumIntensities[spectrumIndex];
				float deltaMass=(float)(targetMass-spectrumMasses[spectrumIndex]);
				peakscores[libraryIndex]=new PeakScores(score, new FragmentIon(targetMass, (byte)0, IonType.y), deltaMass); // FIXME target is a hack!
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
	public float[] auxScore(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return new float[0];
	}
}
