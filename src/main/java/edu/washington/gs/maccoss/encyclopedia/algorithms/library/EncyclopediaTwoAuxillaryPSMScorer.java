package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.algorithms.EncyclopediaAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class EncyclopediaTwoAuxillaryPSMScorer extends EncyclopediaAuxillaryPSMScorer {
	private static final int numPeaksUsedInAverage=3;
	
	private final boolean runXCorr;
	private final LibraryBackgroundInterface background;
	private final SparseXCorrCalculator librarySparseCalculator;
	private final SparseXCorrCalculator sparseModelCalculator;

	public EncyclopediaTwoAuxillaryPSMScorer(SearchParameters parameters, LibraryBackgroundInterface background, boolean runXCorr) {
		super(parameters);
		this.background=background;
		this.runXCorr=runXCorr;
		this.librarySparseCalculator=null;
		this.sparseModelCalculator=null;
	}
	
	
	
	private EncyclopediaTwoAuxillaryPSMScorer(SearchParameters parameters, LibraryBackgroundInterface background, boolean runXCorr, SparseXCorrCalculator librarySparseCalculator, SparseXCorrCalculator sparseModelCalculator) {
		super(parameters);
		this.runXCorr=runXCorr;
		this.background=background;
		this.librarySparseCalculator=librarySparseCalculator;
		this.sparseModelCalculator=sparseModelCalculator;
	}

	@Override
	public EncyclopediaAuxillaryPSMScorer getEntryOptimizedScorer(LibraryEntry entry) {
		SparseXCorrCalculator librarySparse=new SparseXCorrCalculator(entry, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), parameters);
		SparseXCorrCalculator sparseModel=new SparseXCorrCalculator(entry.getPeptideModSeq(), entry.getPrecursorCharge(), parameters);
		return new EncyclopediaTwoAuxillaryPSMScorer(parameters, background, runXCorr, librarySparse, sparseModel);
	}
	
	@Override
	public float[] score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		// precursor scoring
		float[] precursorScores=getPrecursorScores(entry, spectrum.getScanStartTime(), predictedIsotopeDistribution, precursors);
		float averageAbsPPM=precursorScores[0];
		float isotopeDotProduct=precursorScores[1];
		float averagePPM=precursorScores[2];

		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		MassTolerance libraryTolerance=parameters.getLibraryFragmentTolerance();
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		double[] ions=model.getPrimaryIons(parameters.getFragType(), entry.getPrecursorCharge(), false);
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		float[] correlation=entry.getCorrelationArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		float intensityThreshold=spectrum.getTIC()/(1+predictedMasses.length*predictedMasses.length);
		int numberOfMatchingPeaksAboveThreshold=0;
		int numberOfMatchingPeaks=0;
		double dotProduct=0.0;
		double weightedDotProduct=0.0;
		TDoubleArrayList predictedTargets=new TDoubleArrayList();
		TFloatArrayList predictedTargetIntensities=new TFloatArrayList();
		TFloatArrayList actualTargetIntensities=new TFloatArrayList();
		ArrayList<XYPoint> fragmentDeltaMasses=new ArrayList<XYPoint>();
		for (double target : ions) {
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
				if (intensity>0) {
					numberOfMatchingPeaks++;
					if (intensity>intensityThreshold) {
						numberOfMatchingPeaksAboveThreshold++;
					}
				}
				float peakScore=predictedIntensity*intensity*maxCorrelation;
				dotProduct+=peakScore;
				float weight=background==null?1.0f:background.getFraction(target);
				weightedDotProduct+=peakScore*weight;
				predictedTargets.add(target);
				predictedTargetIntensities.add(predictedIntensity);
				actualTargetIntensities.add(intensity);
				
				fragmentDeltaMasses.add(new XYPoint(intensity, deltaMass));
			}
		}
		
		float averageFragmentDeltaMasses=0.0f, averageAbsFragDeltaMass=0.0f;
		if (fragmentDeltaMasses.size()==0) {
			averageAbsFragDeltaMass=(float)acquiredTolerance.getToleranceThreshold();
			averageFragmentDeltaMasses=(float)acquiredTolerance.getToleranceThreshold();
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
				averageAbsFragDeltaMass+=(float)acquiredTolerance.getToleranceThreshold();
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
			float weight=background==null?1.0f:background.getFraction(target);
			weightedSumOfSquaredErrors+=deltaSquared*weight;
		}

		float xTandem;
		if (numberOfMatchingPeaks==0) {
			xTandem=0.0f;
		} else {
			xTandem=((float)Log.protectedLog10(dotProduct))+Log.logFactorial(numberOfMatchingPeaks); // really log10(X!Tandem score)
		}
		
		if (runXCorr) {
			SparseXCorrSpectrum sparseScan=SparseXCorrCalculator.normalize(spectrum, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), false, parameters);
			
			SparseXCorrCalculator librarySparse=librarySparseCalculator!=null?librarySparseCalculator:new SparseXCorrCalculator(entry, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), parameters);
			float xCorrLib=librarySparse.score(sparseScan);
			SparseXCorrCalculator sparseModel=sparseModelCalculator!=null?sparseModelCalculator:new SparseXCorrCalculator(entry.getPeptideModSeq(), entry.getPrecursorCharge(), parameters);
			float xCorrModel=sparseModel.score(sparseScan);
			
			return new float[] {xTandem, xCorrLib, xCorrModel, (float)Log.protectedLog10(dotProduct), (float)Log.protectedLog10(weightedDotProduct), sumOfSquaredErrors, weightedSumOfSquaredErrors, numberOfMatchingPeaks, numberOfMatchingPeaksAboveThreshold, averageAbsFragDeltaMass, averageFragmentDeltaMasses, isotopeDotProduct, averageAbsPPM, averagePPM};
		} else {
			return new float[] {xTandem, (float)Log.protectedLog10(dotProduct), (float)Log.protectedLog10(weightedDotProduct), sumOfSquaredErrors, weightedSumOfSquaredErrors, numberOfMatchingPeaks, numberOfMatchingPeaksAboveThreshold, averageAbsFragDeltaMass, averageFragmentDeltaMasses, isotopeDotProduct, averageAbsPPM, averagePPM};
		}
	}

	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		return getScoreNames(runXCorr);
	}

	public static String[] getScoreNames(boolean runXCorr) {
		// extra scores at the beginning
		if (runXCorr) {
			return new String[] {"e2Score", "evalue", "correlationToGaussian", "correlationToPrecursor", "correlationToPlusOne", "isIntegratedSignal", "isIntegratedPrecursor", "numPeaksWithGoodCorrelation", "numPeaksWithGreatCorrelation", 
					"primary", "xCorrLib", "xCorrModel", "LogDotProduct", "logWeightedDotProduct", "sumOfSquaredErrors", "weightedSumOfSquaredErrors", "numberOfMatchingPeaks", "numberOfMatchingPeaksAboveThreshold", "averageAbsFragmentDeltaMass", "averageFragmentDeltaMasses", "isotopeDotProduct", "averageAbsParentDeltaMass", "averageParentDeltaMass"};
		} else {
			return new String[] {"e2Score", "evalue", "correlationToGaussian", "correlationToPrecursor", "correlationToPlusOne", "isIntegratedSignal", "isIntegratedPrecursor", "numPeaksWithGoodCorrelation", "numPeaksWithGreatCorrelation", 
					"primary", "LogDotProduct", "logWeightedDotProduct", "sumOfSquaredErrors", "weightedSumOfSquaredErrors", "numberOfMatchingPeaks", "numberOfMatchingPeaksAboveThreshold", "averageAbsFragmentDeltaMass", "averageFragmentDeltaMasses", "isotopeDotProduct", "averageAbsParentDeltaMass", "averageParentDeltaMass", "xCorrModel"};
		}
	}
	
	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		float maxFragPPMError=(float)parameters.getFragmentTolerance().getToleranceThreshold();
		float maxPrePPMError=(float)parameters.getPrecursorTolerance().getToleranceThreshold();

		return new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
				maxFragPPMError, maxFragPPMError, 0.0f, maxPrePPMError, maxPrePPMError, 100.0f};
	}
}
