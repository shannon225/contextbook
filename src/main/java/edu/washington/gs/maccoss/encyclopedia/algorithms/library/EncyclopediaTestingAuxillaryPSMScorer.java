package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
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

public class EncyclopediaTestingAuxillaryPSMScorer extends EncyclopediaAuxillaryPSMScorer {
	private static final int numPeaksUsedInAverage=3;
	
	private final LibraryBackgroundInterface background;
	private final SparseXCorrCalculator librarySparseCalculator;
	private final SparseXCorrCalculator sparseModelCalculator;

	public EncyclopediaTestingAuxillaryPSMScorer(SearchParameters parameters, LibraryBackgroundInterface background) {
		super(parameters);
		this.background=background;
		this.librarySparseCalculator=null;
		this.sparseModelCalculator=null;
	}
	
	
	
	private EncyclopediaTestingAuxillaryPSMScorer(SearchParameters parameters, LibraryBackgroundInterface background, SparseXCorrCalculator librarySparseCalculator, SparseXCorrCalculator sparseModelCalculator) {
		super(parameters);
		this.background=background;
		this.librarySparseCalculator=librarySparseCalculator;
		this.sparseModelCalculator=sparseModelCalculator;
	}

	@Override
	public EncyclopediaAuxillaryPSMScorer getEntryOptimizedScorer(LibraryEntry entry) {
		SparseXCorrCalculator librarySparse=new SparseXCorrCalculator(entry, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), parameters);
		SparseXCorrCalculator sparseModel=new SparseXCorrCalculator(entry.getPeptideModSeq(), entry.getPrecursorCharge(), parameters);
		return new EncyclopediaTestingAuxillaryPSMScorer(parameters, background, librarySparse, sparseModel);
	}
	
	@Override
	public float[] score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
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
		
		float[] unnormScores=getScores(predictedTargetIntensitiesArray, actualTargetIntensitiesArray); 
		float[] sumnormScores=getScores(General.normalizeToSum(predictedTargetIntensitiesArray), General.normalizeToSum(actualTargetIntensitiesArray));
		float[] maxnormScores=getScores(General.normalizeToMaxOne(predictedTargetIntensitiesArray), General.normalizeToMaxOne(actualTargetIntensitiesArray));
		float[] l2normScores=getScores(General.normalizeToL2(predictedTargetIntensitiesArray), General.normalizeToL2(actualTargetIntensitiesArray)); 
		
		return General.concatenate(unnormScores, sumnormScores, maxnormScores, l2normScores);
	}

	private float[] getScores(float[] predictedTargetIntensitiesArray, float[] actualTargetIntensitiesArray) {
		float dotProduct=General.sum(General.multiply(predictedTargetIntensitiesArray, actualTargetIntensitiesArray));
		
		float protectedDP=dotProduct;
		if (protectedDP>=1.0f) protectedDP=0.9999f;
		if (Float.isNaN(dotProduct)||dotProduct<=0.0f) return new float[] {0.0f,0.0f,0.0f,0.0f,0.0f};
		
		float contrastAngle=1.0f-(2.0f*(float)Math.acos(protectedDP))/(float)Math.PI;
		float logit=(float)Math.log(dotProduct/(1.0f-protectedDP));
		float sumOfSquaredErrors=0.0f; // normalized to sum of targeted intensities
		
		int numberOfMatchingPeaks=0;
		for (int i=0; i<predictedTargetIntensitiesArray.length; i++) {
			if (predictedTargetIntensitiesArray[i]>0.0&&actualTargetIntensitiesArray[i]>0.0) {
				numberOfMatchingPeaks++;
			
				float delta=predictedTargetIntensitiesArray[i]-actualTargetIntensitiesArray[i];
				float deltaSquared=delta*delta;
				sumOfSquaredErrors+=deltaSquared;
			}
		}

		float xTandem;
		if (dotProduct==0) {
			xTandem=0.0f;
		} else {
			xTandem=((float)Log.protectedLog10(dotProduct))+Log.logFactorial(numberOfMatchingPeaks); // really log10(X!Tandem score)
		}
		
		return new float[] {dotProduct, contrastAngle, logit, -(float)Math.log(sumOfSquaredErrors), xTandem};
	}

	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		//auxScoreArray=General.concatenate(new float[] {scores[1], evalue, correlationToGaussian, 
		//correlationToPrecursor, isIntegratedSignal, isIntegratedPrecursor,
		//numPeaksWithGoodCorrelation}, auxScoreArray);
		return getScoreNames();
	}



	public static String[] getScoreNames() {
		return new String[] {"ScribeScore", "evalue", "correlationToGaussian", "correlationToPrecursor", "isIntegratedSignal", "isIntegratedPrecursor", "numPeaksWithGoodCorrelation", 
				
				"unnormDotProduct", "unnormContrastAngle", "unnormLogit", "unnormScribeScore", "unnormTandem",
				"sumnormDotProduct", "sumnormContrastAngle", "sumnormLogit", "sumnormScribeScore", "sumnormTandem",
				"maxnormDotProduct", "maxnormContrastAngle", "maxnormLogit", "maxnormScribeScore", "maxnormTandem",
				"l2normDotProduct", "l2normContrastAngle", "l2normLogit", "l2normScribeScore", "l2normTandem",};
	}

	@Override
	public int getParentDeltaMassIndex() {
		return AuxillaryPSMScorer.MISSING_INDEX;
	}

	@Override
	public int getFragmentDeltaMassIndex() {
		return AuxillaryPSMScorer.MISSING_INDEX;
	}
	
	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		float maxFragPPMError=(float)parameters.getFragmentTolerance().getToleranceThreshold();
		float maxPrePPMError=(float)parameters.getPrecursorTolerance().getToleranceThreshold();

		return new float[] {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
				maxFragPPMError, maxFragPPMError, 0.0f, maxPrePPMError, maxPrePPMError, 100.0f};
	}
}
