package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class EncyclopediaOneAuxillaryPSMScorer extends AuxillaryPSMScorer {

	public EncyclopediaOneAuxillaryPSMScorer(SearchParameters parameters) {
		super(parameters);
	}
	
	@Override
	public float[] score(LibraryEntry entry, Stripe spectrum, PrecursorScanMap precursors) {
		// precursor scoring
		float[] precursorScores=getPrecursorScores(entry, spectrum.getScanStartTime(), precursors);
		float averageAbsPPM=precursorScores[0];
		float isotopeDotProduct=precursorScores[1];
		float averagePPM=precursorScores[2];

		MassTolerance tolerance=parameters.getFragmentTolerance();
		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		double[] ions=model.getPrimaryIons(parameters.getFragType(), entry.getPrecursorCharge());
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();
		
		int numberOfMatchingPeaks=0;
		float dotProduct=0.0f;
		TFloatArrayList predictedTargetIntensities=new TFloatArrayList();
		TFloatArrayList actualTargetIntensities=new TFloatArrayList();
		TFloatArrayList fragmentDeltaMasses=new TFloatArrayList();
		float averageAbsFragDeltaMass=0.0f;
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
						deltaMass=(float)((target-predictedMasses[indicies[j]])*1000000.0/target);
					}
				}
				if (intensity>0) {
					numberOfMatchingPeaks++;
				}
				dotProduct+=predictedIntensity*intensity;
				predictedTargetIntensities.add(predictedIntensity);
				actualTargetIntensities.add(intensity);
				fragmentDeltaMasses.add(deltaMass);
				averageAbsFragDeltaMass+=Math.abs(deltaMass);
			}
		}
		averageAbsFragDeltaMass=averageAbsFragDeltaMass/fragmentDeltaMasses.size();
		
		float averageFragmentDeltaMasses=General.mean(fragmentDeltaMasses.toArray());

		float[] predictedTargetIntensitiesArray=predictedTargetIntensities.toArray();
		float[] actualTargetIntensitiesArray=actualTargetIntensities.toArray();
		
		float sumPredictedTargets=General.sum(predictedTargetIntensitiesArray);
		float sumActualTargets=General.sum(actualTargetIntensitiesArray);
		
		float sumOfSquaredErrors=0.0f; // normalized to sum of targeted intensities
		for (int i=0; i<predictedTargetIntensitiesArray.length; i++) {
			float predicted=predictedTargetIntensitiesArray[i]/sumPredictedTargets;
			float actual=actualTargetIntensitiesArray[i]/sumActualTargets;
			float delta=predicted-actual;
			sumOfSquaredErrors+=delta*delta;
		}
		
		return new float[] {dotProduct, sumOfSquaredErrors, numberOfMatchingPeaks, averageAbsFragDeltaMass, averageFragmentDeltaMasses, isotopeDotProduct, averageAbsPPM, averagePPM};
	}

	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		return new String[] {"dotProduct", "sumOfSquaredErrors", "numberOfMatchingPeaks", "averageAbsFragDeltaMass", "averageFragmentDeltaMasses", "isotopeDotProduct", "averageAbsPPM", "averagePPM"};
	}
	
	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		float maxFragPPMError=(float)parameters.getFragmentTolerance().getPpmTolerance();
		float maxPrePPMError=(float)parameters.getPrecursorTolerance().getPpmTolerance();
	
		return new float[] {0.0f, 0.0f, maxFragPPMError, maxFragPPMError, 0.0f, maxPrePPMError, maxPrePPMError};
	}
}
