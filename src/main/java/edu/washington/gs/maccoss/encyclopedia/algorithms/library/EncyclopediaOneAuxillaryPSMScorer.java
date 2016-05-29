package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class EncyclopediaOneAuxillaryPSMScorer extends AuxillaryPSMScorer {
	private static final int numPeaksUsedInAverage=3;
	
	private final LibraryBackground background;

	public EncyclopediaOneAuxillaryPSMScorer(SearchParameters parameters, LibraryBackground background) {
		super(parameters);
		this.background=background;
	}
	
	@Override
	public float[] score(LibraryEntry entry, Stripe spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		// precursor scoring
		float[] precursorScores=getPrecursorScores(entry, spectrum.getScanStartTime(), predictedIsotopeDistribution, precursors);
		float averageAbsPPM=precursorScores[0];
		float isotopeDotProduct=precursorScores[1];
		float averagePPM=precursorScores[2];

		MassTolerance tolerance=parameters.getFragmentTolerance();
		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		double[] ions=model.getPrimaryIons(parameters.getFragType(), entry.getPrecursorCharge());
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		float[] correlation=entry.getCorrelationArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();
		
		int numberOfMatchingPeaks=0;
		double dotProduct=0.0;
		double weightedDotProduct=0.0;
		TDoubleArrayList predictedTargets=new TDoubleArrayList();
		TFloatArrayList predictedTargetIntensities=new TFloatArrayList();
		TFloatArrayList actualTargetIntensities=new TFloatArrayList();
		ArrayList<XYPoint> fragmentDeltaMasses=new ArrayList<XYPoint>();
		for (double target : ions) {
			int[] predictedIndicies=tolerance.getIndicies(predictedMasses, target);
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
				if (intensity>0) {
					numberOfMatchingPeaks++;
				}
				float peakScore=predictedIntensity*intensity*maxCorrelation;
				dotProduct+=peakScore;
				float fraction=background.getFraction(target);
				weightedDotProduct+=peakScore*fraction;
				predictedTargets.add(target);
				predictedTargetIntensities.add(predictedIntensity);
				actualTargetIntensities.add(intensity);
				
				fragmentDeltaMasses.add(new XYPoint(intensity, deltaMass));
			}
		}
		
		float averageFragmentDeltaMasses=0.0f, averageAbsFragDeltaMass=0.0f;
		if (fragmentDeltaMasses.size()==0) {
			averageAbsFragDeltaMass=(float)tolerance.getPpmTolerance();
			averageFragmentDeltaMasses=(float)tolerance.getPpmTolerance();
		} else {
			Collections.sort(fragmentDeltaMasses);
			Collections.reverse(fragmentDeltaMasses);
			
			int count=0;
			for (XYPoint xyPoint : fragmentDeltaMasses) {
				averageFragmentDeltaMasses+=(float)xyPoint.y;
				averageAbsFragDeltaMass+=Math.abs((float)xyPoint.y);
				count++;
				if (count>numPeaksUsedInAverage) break;
			}
			for (int i=count; i<numPeaksUsedInAverage; i++) {
				averageAbsFragDeltaMass+=(float)tolerance.getPpmTolerance();
			}
			averageFragmentDeltaMasses=averageFragmentDeltaMasses/count;
			averageAbsFragDeltaMass=averageAbsFragDeltaMass/numPeaksUsedInAverage;
		}

		float[] predictedTargetIntensitiesArray=predictedTargetIntensities.toArray();
		float[] actualTargetIntensitiesArray=actualTargetIntensities.toArray();
		
		float sumPredictedTargets=General.sum(predictedTargetIntensitiesArray);
		float sumActualTargets=General.sum(actualTargetIntensitiesArray);
		
		float sumOfSquaredErrors=0.0f; // normalized to sum of targeted intensities
		float weightedSumOfSquaredErrors=0.0f;

		if (predictedTargetIntensitiesArray.length==0) {
			sumOfSquaredErrors=1.0f;
			weightedSumOfSquaredErrors=10.0f;
		}
		
		for (int i=0; i<predictedTargetIntensitiesArray.length; i++) {
			float predicted=predictedTargetIntensitiesArray[i]/sumPredictedTargets;
			float actual;
			if (sumActualTargets==0.0f) {
				actual=0.0f;
			} else {
				actual=actualTargetIntensitiesArray[i]/sumActualTargets;
			}
			float delta=predicted-actual;
			float deltaSquared=delta*delta;
			double target=predictedTargets.get(i);
			sumOfSquaredErrors+=deltaSquared;
			weightedSumOfSquaredErrors+=deltaSquared*background.getFraction(target);
		}

		float xTandem;
		if (numberOfMatchingPeaks==0) {
			xTandem=0.0f;
		} else {
			xTandem=((float)Log.protectedLog10(dotProduct))+Log.logFactorial(numberOfMatchingPeaks); // really log10(X!Tandem score)
		}
		
		return new float[] {xTandem, (float)Log.protectedLog10(dotProduct), (float)Log.protectedLog10(weightedDotProduct), sumOfSquaredErrors, weightedSumOfSquaredErrors, numberOfMatchingPeaks, averageAbsFragDeltaMass, averageFragmentDeltaMasses, isotopeDotProduct, averageAbsPPM, averagePPM};
	}

	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		return getScoreNames();
	}

	public static String[] getScoreNames() {
		return new String[] {"xTandem", "LogDotProduct", "logWeightedDotProduct", "sumOfSquaredErrors", "weightedSumOfSquaredErrors", "numberOfMatchingPeaks", "averageAbsFragDeltaMass", "averageFragmentDeltaMasses", "isotopeDotProduct", "averageAbsPPM", "averagePPM", "eValue"};
	}
	
	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		float maxFragPPMError=(float)parameters.getFragmentTolerance().getPpmTolerance();
		float maxPrePPMError=(float)parameters.getPrecursorTolerance().getPpmTolerance();
	
		return new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, maxFragPPMError, maxFragPPMError, 0.0f, maxPrePPMError, maxPrePPMError};
	}
}
