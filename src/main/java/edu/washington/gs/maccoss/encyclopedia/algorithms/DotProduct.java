package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AnnotatedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
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
		if (entry instanceof AnnotatedLibraryEntry) {
			return getIndividualPeakScores(entry, spectrum, normalize, ((AnnotatedLibraryEntry)entry).getIonAnnotations());
		} else {
			throw new EncyclopediaException("DotProduct doesn't currently handle ion selection. Please report this bug!");
		}
	}
	
	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize, FragmentIon[] ions) {
		MassTolerance acquiredTolerance=tolerance;
		MassTolerance libraryTolerance=tolerance;
		
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
						deltaMass=(float)((target-acquiredMasses[indicies[j]])*1000000.0/target);
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
	public float[] auxScore(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return new float[0];
	}
}
