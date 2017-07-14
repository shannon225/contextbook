package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.*;

import java.util.ArrayList;

//@Immutable
public class DotProduct implements PSMPeakScorer {
	private final SearchParameters parameters;

	public DotProduct(SearchParameters parameters) {
		this.parameters=parameters;
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer#score(edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry, edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe)
	 */
	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return PeakScores.sumScores(getIndividualPeakScores(entry, spectrum, false));
	}
	
	@Override
	public float score(LibraryEntry entry, Spectrum spectrum, FragmentIon[] ions) {
		return PeakScores.sumScores(getIndividualPeakScores(entry, spectrum, false));
	}
	
	@Override
	public String[] getAuxScoreNames(LibraryEntry entry) {
		return new String[0];
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize) {
		if (entry instanceof AnnotatedLibraryEntry) {
			return getIndividualPeakScores(entry, spectrum, normalize, ((AnnotatedLibraryEntry)entry).getIonAnnotations());
		} else {
			double[] mzs=entry.getMassArray();
			FragmentIon[] ions=new FragmentIon[mzs.length];
			for (byte i=0; i<ions.length; i++) {
				ions[i]=new FragmentIon(mzs[i], (byte)(i+1), IonType.unknown);
			}
			
			return getIndividualPeakScores(entry, spectrum, normalize, ions);
		}
	}
	
	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Spectrum spectrum, boolean normalize, FragmentIon[] ions) {
		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		MassTolerance libraryTolerance=parameters.getLibraryFragmentTolerance();
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		float[] correlation=entry.getCorrelationArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		ArrayList<PeakScores> scoredPeaks=new ArrayList<PeakScores>();
		for (FragmentIon targetIon : ions) {
			double target=targetIon.mass;
			
			int[] predictedIndicies=libraryTolerance.getIndicies(predictedMasses, target);
			float predictedIntensity=0.0f;
			float maxCorrelation=0.01f;
			for (int i=0; i<predictedIndicies.length; i++) {
				if (predictedIntensity<predictedIntensities[predictedIndicies[i]]) {
					predictedIntensity=predictedIntensities[predictedIndicies[i]];
				}
				if (maxCorrelation<correlation[predictedIndicies[i]]) {
					maxCorrelation=correlation[predictedIndicies[i]];
				}
			}
			
			if (predictedIntensity>0) {
				int[] indicies=acquiredTolerance.getIndicies(acquiredMasses, target);
				float intensity=0.0f;
				float bestPeakIntensity=0.0f;
				float deltaMass=0.0f;
				for (int j=0; j<indicies.length; j++) {
					intensity+=acquiredIntensities[indicies[j]];
					
					if (acquiredIntensities[indicies[j]]>bestPeakIntensity) {
						bestPeakIntensity=acquiredIntensities[indicies[j]];

						deltaMass=(float)acquiredTolerance.getDeltaScore(target, acquiredMasses[indicies[j]]);
					}
				}
				float peakScore=predictedIntensity*intensity*maxCorrelation;
				if (intensity>0.0f) {
					scoredPeaks.add(new PeakScores(peakScore, targetIon, deltaMass));
				} else {
					scoredPeaks.add(null);
				}
			}
		}
		return scoredPeaks.toArray(new PeakScores[scoredPeaks.size()]);
	}
	
	@Override
	public float[] auxScore(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return new float[0];
	}
}
