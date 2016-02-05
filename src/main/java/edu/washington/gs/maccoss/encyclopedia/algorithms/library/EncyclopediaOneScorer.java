package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;

public class EncyclopediaOneScorer implements PSMScorer {
	private final SearchParameters parameters;
	private final EncyclopediaOneAuxillaryPSMScorer auxScorer;

	public EncyclopediaOneScorer(SearchParameters parameters, LibraryBackground background) {
		this.parameters=parameters;
		auxScorer=new EncyclopediaOneAuxillaryPSMScorer(parameters, background);
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
	public float score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		return score(entry, spectrum);
	}

	public float score(LibraryEntry entry, Stripe spectrum) {
		PeakScores[] individualPeakScores=getIndividualPeakScores(entry, spectrum, true);
		int count=0; // number of matches
		for (int i=0; i<individualPeakScores.length; i++) {
			if (individualPeakScores[i]!=null) count++;
		}
		
		if (count==0) return 0.0f;
		
		float dotProduct=PeakScores.sumScores(individualPeakScores); // dot product
		
		return Log.protectedLog10(dotProduct)+Log.logFactorial(count); // X!Tandem score
	}

	@Override
	public PeakScores[] getIndividualPeakScores(LibraryEntry entry, Stripe spectrum, boolean normalize) {
		MassTolerance tolerance=parameters.getFragmentTolerance();
		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		double[] ions=model.getPrimaryIons(parameters.getFragType(), entry.getPrecursorCharge());
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();
		
		ArrayList<PeakScores> scoredPeaks=new ArrayList<PeakScores>();
		for (double target : ions) {
			float predictedIntensity=tolerance.getIntegratedIntensity(predictedMasses, predictedIntensities, target);
			
			if (predictedIntensity>0) {
				int[] indicies=tolerance.getIndicies(acquiredMasses, target);
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
				float peakScore=predictedIntensity*intensity;
				if (intensity>0.0f) {
					scoredPeaks.add(new PeakScores(peakScore, target, deltaMass));
				} else {
					scoredPeaks.add(null);
				}
			}
		}
		return scoredPeaks.toArray(new PeakScores[scoredPeaks.size()]);
	}

}
