package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
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

public class ScribeAuxillaryPSMScorer extends AuxillaryPSMScorer {
	private static final int numPeaksUsedInAverage=3;

	public ScribeAuxillaryPSMScorer(SearchParameters parameters) {
		super(parameters);
	}
	
	@Override
	public float[] score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		// precursor scoring
		float[] precursorScores=getPrecursorScores(entry, spectrum.getScanStartTime(), predictedIsotopeDistribution, precursors);
		float averageAbsPPM=precursorScores[0];
		float isotopeDotProduct=precursorScores[1];
		float averagePPM=precursorScores[2];
		float percentBlankOverMono=precursorScores[3];
		float numberPrecursorMatch=precursorScores[4];

		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		MassTolerance libraryTolerance=parameters.getLibraryFragmentTolerance();
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		double[] ions=model.getPrimaryIons(parameters.getFragType(), entry.getPrecursorCharge(), false);
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		float[] correlation=entry.getCorrelationArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		int numberOfMatchingPeaks=0;
		double dotProduct=0.0;
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
				}
				float peakScore=predictedIntensity*intensity*maxCorrelation;
				dotProduct+=peakScore;
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

		if (predictedTargetIntensitiesArray.length==0) {
			sumOfSquaredErrors=1.0f;
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
			sumOfSquaredErrors+=deltaSquared;
		}

		float xTandem;
		if (numberOfMatchingPeaks==0) {
			xTandem=0.0f;
		} else {
			xTandem=((float)Log.protectedLog10(dotProduct))+Log.logFactorial(numberOfMatchingPeaks); // really log10(X!Tandem score)
		}
		
		SparseXCorrSpectrum sparseScan=SparseXCorrCalculator.normalize(spectrum, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), false, parameters);
		
		SparseXCorrCalculator librarySparse=new SparseXCorrCalculator(entry, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), parameters);
		float xCorrLib=librarySparse.score(sparseScan);
		SparseXCorrCalculator sparseModel=new SparseXCorrCalculator(entry.getPeptideModSeq(), entry.getPrecursorCharge(), parameters);
		float xCorrModel=sparseModel.score(sparseScan);
			
		return new float[] {(float)Log.protectedLog10(dotProduct), xTandem, xCorrLib, xCorrModel, Log.protectedLn(1.0f/sumOfSquaredErrors), 
				numberOfMatchingPeaks, averageAbsFragDeltaMass, averageFragmentDeltaMasses, isotopeDotProduct, 
				averageAbsPPM, averagePPM, percentBlankOverMono, numberPrecursorMatch};
	}

	public static String[] getScoreNames() {
		return new String[] {"LogDotProduct", "primary", "xCorrLib", "xCorrModel", "lnInvSumOfSquaredErrors", 
				"numberOfMatchingPeaks", "averageAbsFragmentDeltaMass", "averageFragmentDeltaMasses", "isotopeDotProduct", 
				"averageAbsParentDeltaMass", "averageParentDeltaMass", "percentBlankOverMono", "numberPrecursorMatch", "eValue", "numConsidered"};
	}
	
	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		return getScoreNames();
	}
	
	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		throw new EncyclopediaException("Sorry, missing score data is not allowed for Scribe. Please contact the help boards if you see this message.");
	}
}
