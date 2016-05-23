package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IonType;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;

//@Immutable
public class PecanRawScorer implements PSMScorer {
	private final MassTolerance fragmentTolerance;
	private final AuxillaryPSMScorer auxScorer;

	public PecanRawScorer(MassTolerance fragmentTolerance, AuxillaryPSMScorer auxScorer) {
		this.fragmentTolerance = fragmentTolerance;
		this.auxScorer=auxScorer;
	}

	@Override
	public float score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return score(entry, spectrum);
	}

	public float score(LibraryEntry entry, Stripe spectrum) {
		return PeakScores.sumScores(getIndividualPeakScores(entry, spectrum, true)); // dot product
	}

	@Override
	public float[] auxScore(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return auxScorer.score(entry, spectrum, predictedIsotopeDistribution, precursors);
	}
	
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return auxScorer.getScoreNames(entry);
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize) {
		return getIndividualPeakScores(entry, spectrum, normalize, null);
	}
	
	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize, FragmentIon[] ions) {
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
			int[] indicies=fragmentTolerance.getIndicies(spectrumMasses, targetMass);
			float intensity=0.0f;
			float bestPeakIntensity=0.0f;
			float deltaMass=0.0f;
			for (int j=0; j<indicies.length; j++) {
				intensity+=spectrumIntensities[indicies[j]];
				
				if (spectrumIntensities[indicies[j]]>bestPeakIntensity) {
					bestPeakIntensity=spectrumIntensities[indicies[j]];
					deltaMass=(float)((targetMass-spectrumMasses[indicies[j]])*1000000.0/targetMass);
				}
			}
			float peakScore=libraryIntensities[i]*intensity;
			if (intensity>0.0f) {
				individualPeakScores[i]=new PeakScores(peakScore, new FragmentIon(targetMass, (byte)0, IonType.y), deltaMass);// FIXME target is a hack!
			}
		}
		
		return individualPeakScores;
	}
}
