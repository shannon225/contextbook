package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

//@Immutable
public class PecanRawScorer implements PSMPeakScorer {
	public static final String PRIMARY_SCORE_NAME="peakCalibratedScore";
	
	private final MassTolerance fragmentTolerance;
	private final AuxillaryPSMScorer auxScorer;

	public PecanRawScorer(MassTolerance fragmentTolerance, AuxillaryPSMScorer auxScorer) {
		this.fragmentTolerance = fragmentTolerance;
		this.auxScorer=auxScorer;
	}
	
	public static String getPrimaryScoreName() {
		return PRIMARY_SCORE_NAME;
	}

	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return score(entry, spectrum);
	}

	public float score(LibraryEntry entry, Spectrum spectrum) {
		return PeakScores.sumScores(getIndividualPeakScores(entry, spectrum, true)); // dot product
	}
	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, FragmentIon[] ions) {
		return PeakScores.sumScores(getIndividualPeakScores(entry, spectrum, true)); // dot product
	}

	@Override
	public float[] auxScore(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return auxScorer.score(entry, spectrum, predictedIsotopeDistribution, precursors);
	}
	
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return auxScorer.getScoreNames(entry);
	}

	@Override
	public int getParentDeltaMassIndex() {
		return auxScorer.getParentDeltaMassIndex();
	}

	@Override
	public int getFragmentDeltaMassIndex() {
		return auxScorer.getFragmentDeltaMassIndex();
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize) {
		//if (entry instanceof AnnotatedLibraryEntry) {
		//	return getIndividualPeakScores(entry, spectrum, normalize, ((AnnotatedLibraryEntry) entry).getIonAnnotations());
		//} else {
			return getIndividualPeakScores(entry, spectrum, normalize, null);
		//}
	}
	
	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize, FragmentIon[] ions) {
		if (ions!=null) {
			throw new EncyclopediaException("PecanRawScorer doesn't currently handle ion selection. Please report this bug!");
		}
		
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
		
		PeakScores[] individualPeakScores=new PeakScores[libraryMasses.length];
		
		if (libraryMasses.length==0||spectrumMasses.length==0) return individualPeakScores;

		for (int i=0; i<libraryMasses.length; i++) {
			double targetMass=libraryMasses[i];
			int[] indicies=fragmentTolerance.getIndices(spectrumMasses, targetMass);
			float intensity=0.0f;
			float bestPeakIntensity=0.0f;
			float deltaMass=0.0f;
			for (int j=0; j<indicies.length; j++) {
				// NOTE, this is the sum of square rooted intensities, as done by python-pecan
				intensity+=spectrumIntensities[indicies[j]];
				
				if (spectrumIntensities[indicies[j]]>bestPeakIntensity) {
					bestPeakIntensity=spectrumIntensities[indicies[j]];

					// NOTE, the error is derived from the apex, which is NOT
					// done by python-pecan. Python-pecan currently takes the
					// minimum delta mass of all peaks in the window, which
					// Sonia agrees isn't correct.
					deltaMass=(float)fragmentTolerance.getDeltaScore(targetMass, spectrumMasses[indicies[j]]);
				}
			}
			float peakScore=libraryIntensities[i]*intensity;
			if (intensity>0.0f) {
				individualPeakScores[i]=new PeakScores(peakScore, new FragmentIon(targetMass, (byte)0, IonType.y), deltaMass);// TODO target ion index is a hack!
			}
		}
		
		return individualPeakScores;
	}
}
