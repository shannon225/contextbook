package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;
import gnu.trove.procedure.TCharDoubleProcedure;

public class ScribeAuxillaryPSMScorer extends AuxillaryPSMScorer {
	private static final int numPeaksUsedInAverage=3;
	private static final TCharDoubleHashMap targetImmoniumIons=new TCharDoubleHashMap();
	{
		targetImmoniumIons.put('H', 110.0718);
		targetImmoniumIons.put('Y', 136.0762);
		targetImmoniumIons.put('W', 159.0922);
		targetImmoniumIons.put('M', 104.0534);
		targetImmoniumIons.put('F', 120.0813);
	}

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
		FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false);
		HashMap<IonType, float[]> matchedIonLadders=new HashMap<>();
		
		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		float[] correlation=entry.getCorrelationArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		int numberOfMatchingPeaks=0;
		TDoubleArrayList predictedTargets=new TDoubleArrayList();
		TFloatArrayList predictedTargetIntensities=new TFloatArrayList();
		TFloatArrayList actualTargetIntensities=new TFloatArrayList();
		ArrayList<XYPoint> fragmentDeltaMasses=new ArrayList<XYPoint>();
		for (FragmentIon target : ions) {
			double targetMass=target.getMass();
			int[] predictedIndicies=libraryTolerance.getIndicies(predictedMasses, targetMass);
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
				int[] indicies=acquiredTolerance.getIndicies(acquiredMasses, targetMass);
				float intensity=0.0f;
				float bestPeakIntensity=0.0f;
				float deltaMass=0.0f;
				for (int j=0; j<indicies.length; j++) {
					intensity+=acquiredIntensities[indicies[j]];
					
					if (acquiredIntensities[indicies[j]]>bestPeakIntensity) {
						bestPeakIntensity=acquiredIntensities[indicies[j]];

						deltaMass=(float)acquiredTolerance.getDeltaScore(targetMass, acquiredMasses[indicies[j]]);
					}
				}
				if (intensity>0) {
					numberOfMatchingPeaks++;
					
					float[] ladder=matchedIonLadders.get(target.getType());
					if (ladder==null) {
						ladder=new float[entry.getPeptideSeq().length()];
						matchedIonLadders.put(target.getType(), ladder);
					}
					if (target.getIndex()<=ladder.length) {
						ladder[target.getIndex()-1]=intensity;
					}
				}
				predictedTargets.add(targetMass);
				predictedTargetIntensities.add(predictedIntensity);
				actualTargetIntensities.add(intensity);
				
				fragmentDeltaMasses.add(new XYPoint(intensity, deltaMass));
			}
		}
		
		int maxLadderLength=0;
		for (float[] ladder : matchedIonLadders.values()) {
			int currentLength=0;
			for (int i = 0; i < ladder.length; i++) {
				if (ladder[i]>0.0f) {
					currentLength++;
				}
				if (currentLength>maxLadderLength) {
					maxLadderLength=currentLength;
				}
			}
		}
		
		int[] numImmoniumIonsFound=new int[1];
		targetImmoniumIons.forEachEntry(new TCharDoubleProcedure() {
			
			@Override
			public boolean execute(char aa, double targetMass) {
				if (entry.getPeptideSeq().indexOf(aa)>=0) {
					int[] indicies=acquiredTolerance.getIndicies(acquiredMasses, targetMass);

					for (int j=0; j<indicies.length; j++) {
						if (acquiredIntensities[indicies[j]]>0) {
							numImmoniumIonsFound[0]++;
							break;
						};
					}
				}
				return true;
			}
		});
		
		float beta=0.075f*maxLadderLength;
		float rho=0.15f*numImmoniumIonsFound[0];
		float sp=General.sum(actualTargetIntensities.toArray())*numberOfMatchingPeaks*(1.0f+beta)*(1.0f+rho)/ions.length;
		
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

		float[] predictedTargetIntensitiesArray=General.normalizeToL2(predictedTargetIntensities.toArray());
		float[] actualTargetIntensitiesArray=General.normalizeToL2(actualTargetIntensities.toArray());
		
		float dotProduct=General.sum(General.multiply(predictedTargetIntensitiesArray, actualTargetIntensitiesArray));

		if (Float.isNaN(dotProduct)||dotProduct<0.0f) dotProduct=0.0f;
		float protectedDP=dotProduct;
		if (protectedDP>=1.0f) protectedDP=0.99999f;
		if (protectedDP<=0.0f) protectedDP=0.00001f;
		
		float contrastAngle=1.0f-(2.0f*(float)Math.acos(protectedDP))/(float)Math.PI;
		float logit=(float)Math.log(protectedDP/(1.0f-protectedDP));
		float sumOfSquaredErrors=0.0f; // normalized to sum of targeted intensities
		
		for (int i=0; i<predictedTargetIntensitiesArray.length; i++) {
			if (predictedTargetIntensitiesArray[i]>0.0&&actualTargetIntensitiesArray[i]>0.0) {
				float delta=predictedTargetIntensitiesArray[i]-actualTargetIntensitiesArray[i];
				float deltaSquared=delta*delta;
				sumOfSquaredErrors+=deltaSquared;
			}
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
			
		return new float[] {dotProduct, contrastAngle, logit, xTandem, xCorrLib, xCorrModel, Log.protectedLn(1.0f/sumOfSquaredErrors), 
				numberOfMatchingPeaks, averageAbsFragDeltaMass, averageFragmentDeltaMasses, isotopeDotProduct, 
				averageAbsPPM, averagePPM, percentBlankOverMono, numberPrecursorMatch, sp, maxLadderLength};
	}

	public static String[] getScoreNames() {
		return new String[] {"DotProduct", "contrastAngle", "logit", "primary", "xCorrLib", "xCorrModel", "lnInvSumOfSquaredErrors", 
				"numberOfMatchingPeaks", "averageAbsFragmentDeltaMass", "averageFragmentDeltaMasses", "isotopeDotProduct", 
				"averageAbsParentDeltaMass", "averageParentDeltaMass", "percentBlankOverMono", "numberPrecursorMatch", "Sp", "maxLadderLength", 
				"eValue", "numConsidered", "deltaCn", "chargeMatch"};
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
