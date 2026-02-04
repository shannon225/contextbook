package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.EncyclopediaAuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
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

public class EncyclopediaTwoPointOneAuxillaryPSMScorer extends EncyclopediaAuxillaryPSMScorer {
	private static final int numPeaksUsedInAverage=3;
	private final float maxPPMError;
	
	private static final TCharDoubleHashMap targetImmoniumIons=new TCharDoubleHashMap();
	{
		targetImmoniumIons.put('H', 110.0718);
		targetImmoniumIons.put('Y', 136.0762);
		targetImmoniumIons.put('W', 159.0922);
		targetImmoniumIons.put('M', 104.0534);
		targetImmoniumIons.put('F', 120.0813);
	}
	
	private final SparseXCorrCalculator librarySparseCalculator;
	private final SparseXCorrCalculator sparseModelCalculator;

	public EncyclopediaTwoPointOneAuxillaryPSMScorer(SearchParameters parameters) {
		super(parameters);
		this.librarySparseCalculator=null;
		this.sparseModelCalculator=null;
		maxPPMError=(float)parameters.getPrecursorTolerance().getToleranceThreshold();
	}
	
	private EncyclopediaTwoPointOneAuxillaryPSMScorer(SearchParameters parameters, SparseXCorrCalculator librarySparseCalculator, SparseXCorrCalculator sparseModelCalculator) {
		super(parameters);
		this.librarySparseCalculator=librarySparseCalculator;
		this.sparseModelCalculator=sparseModelCalculator;
		maxPPMError=(float)parameters.getPrecursorTolerance().getToleranceThreshold();
	}

	@Override
	public EncyclopediaAuxillaryPSMScorer getEntryOptimizedScorer(LibraryEntry entry) {
		SparseXCorrCalculator librarySparse=new SparseXCorrCalculator(entry, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), parameters);
		SparseXCorrCalculator sparseModel=new SparseXCorrCalculator(entry.getPeptideModSeq(), entry.getPrecursorCharge(), parameters);
		return new EncyclopediaTwoPointOneAuxillaryPSMScorer(parameters, librarySparse, sparseModel);
	}
	
	@Override
	public float[] score(LibraryEntry entry, Spectrum spectrum, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		// precursor scoring
		float[] precursorScores=getLocalPrecursorScores(entry, spectrum.getScanStartTime(), predictedIsotopeDistribution, precursors);
		float isotopeDotProduct=precursorScores[0];
		float averagePPM=precursorScores[1];

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

		float intensityThreshold=spectrum.getTIC()/(1+predictedMasses.length*predictedMasses.length);
		int numberOfMatchingPeaksAboveThreshold=0;
		int numberOfMatchingPeaks=0;
		TDoubleArrayList predictedTargets=new TDoubleArrayList();
		TFloatArrayList predictedTargetIntensities=new TFloatArrayList();
		TFloatArrayList actualTargetIntensities=new TFloatArrayList();
		ArrayList<XYPoint> fragmentDeltaMasses=new ArrayList<XYPoint>();
		for (FragmentIon target : ions) {
			double targetMass=target.getMass();
			int[] predictedIndicies=libraryTolerance.getIndices(predictedMasses, targetMass);
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
				int[] indicies=acquiredTolerance.getIndices(acquiredMasses, targetMass);
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
					if (intensity>intensityThreshold) {
						numberOfMatchingPeaksAboveThreshold++;
					}
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
					int[] indicies=acquiredTolerance.getIndices(acquiredMasses, targetMass);

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
		if (sp<1) sp=1;
		
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
		
		float sumOfSquaredErrors=0.0f; // normalized to sum of targeted intensities

		if (predictedTargetIntensitiesArray.length==0) {
			sumOfSquaredErrors=1.0f;
		} else {
			for (int i=0; i<predictedTargetIntensitiesArray.length; i++) {
				if (predictedTargetIntensitiesArray[i]>0.0&&actualTargetIntensitiesArray[i]>0.0) {
					float delta=predictedTargetIntensitiesArray[i]-actualTargetIntensitiesArray[i];
					float deltaSquared=delta*delta;
					sumOfSquaredErrors+=deltaSquared;
				}
			}
		}
		
		float scribe;
		if (sumOfSquaredErrors<=0.0f) {
			scribe=0.0f;
		} else {
			// shifted +1 to capture part of the negative range while ensuring a non-negative score
			scribe=Log.protectedLn(1.0f/sumOfSquaredErrors)+1f;
		}
		if (scribe < 0.0f) {
			scribe = 0.0f;
		}
		
		float xTandem;
		if (numberOfMatchingPeaks==0||dotProduct<=0) {
			xTandem=0.0f;
		} else {
			// really log10(X!Tandem score)
			// shifted +2 to capture part of the negative range while ensuring a non-negative score
			xTandem=((float)Log.protectedLog10(dotProduct))+Log.logFactorial(numberOfMatchingPeaks)+2f;
		}

		if (xTandem < 0.0f) {
			xTandem = 0.0f;
		}
		
		SparseXCorrSpectrum sparseScan=SparseXCorrCalculator.normalize(spectrum, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), false, parameters);
		
		SparseXCorrCalculator librarySparse=librarySparseCalculator!=null?librarySparseCalculator:new SparseXCorrCalculator(entry, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), parameters);
		float xCorrLib=librarySparse.score(sparseScan);
		SparseXCorrCalculator sparseModel=sparseModelCalculator!=null?sparseModelCalculator:new SparseXCorrCalculator(entry.getPeptideModSeq(), entry.getPrecursorCharge(), parameters);
		float xCorrModel=sparseModel.score(sparseScan);
		
		return new float[] {xTandem, xCorrLib, xCorrModel, scribe, numberOfMatchingPeaks, numberOfMatchingPeaksAboveThreshold, 
				averageFragmentDeltaMasses, isotopeDotProduct, averagePPM, Log.protectedLn(sp), 
				maxLadderLength};
	
	}

	@Override
	public String[] getScoreNames(LibraryEntry entry) {
		return getScoreNames();
	}

	public static String[] getScoreNames() {
		return new String[] {"primary", "evalue", "correlationToGaussian", "correlationToPrecursor", "isIntegratedSignal", "isIntegratedPrecursor", 
				"numPeaksWithGoodCorrelation", "HyperScore", "xCorrLib", "xCorrModel", "scribe", "numberOfMatchingPeaks", 
				"numberOfMatchingPeaksAboveThreshold", "averageFragmentDeltaMasses", "isotopeDotProduct", 
				"averageParentDeltaMass", "lnSp", "maxLadderLength"};
	}
	
	@Override
	public int getFragmentDeltaMassIndex() {
		return 6; // in score array, not in names
	}
	
	@Override
	public int getParentDeltaMassIndex() {
		return 8; // in score array, not in names
	}
	
	@Override
	public float[] getMissingDataScores(LibraryEntry entry) {
		throw new EncyclopediaException("Sorry, missing score data is not allowed for EncyclopeDIA v2 Scoring. Please contact the help boards if you see this message.");
	}
	
	public float[] getLocalPrecursorScores(LibraryEntry entry, float spectrumRT, float[] predictedIsotopeDistribution, PrecursorScanMap precursors) {
		byte charge=entry.getPrecursorCharge();
		Peak[] precursorPacket=precursors.getIsotopePacket(entry.getPrecursorMZ(), spectrumRT, charge, parameters.getPrecursorTolerance());
		if (precursorPacket.length==0) {
			return new float[] {0.0f, maxPPMError};
		}
		
		Triplet<double[], float[], Optional<float[]>> pair=Peak.toArrays(precursorPacket);
		double[] masses=pair.x;
		float[] intensities=pair.y;
		
		// weighted average AbsPPM
		float averagePPM=0.0f; // FINAL SCORE
		// start at 1 to drop "-1" isotope
		float denominator=0.0f;
		for (int i = 1; i < masses.length; i++) {
			byte isotope=(byte)(i-1);
			double predictedMass=entry.getPrecursorMZ()+(isotope*MassConstants.neutronMass/charge);
			float predictedIntensity=predictedIsotopeDistribution[isotope];
			
			if (intensities[i]>0) {
				averagePPM+=(float)parameters.getPrecursorTolerance().getDeltaScore(predictedMass, masses[i])*predictedIntensity;
				denominator+=predictedIntensity;
			}
		}
		
		if (denominator>0) {
			averagePPM=averagePPM/denominator;
		} else {
			return new float[] {0.0f, maxPPMError};
		}
		
		// precursor idotp
		intensities=IsotopicDistributionCalculator.normalizeToMax(intensities);
		float isotopeDotProduct=0.0f; // FINAL SCORE
		float euclideanDistanceIntensities=0.0f;
		float euclideanDistancePredicted=0.0f;
		for (int i = 0; i < PrecursorScanMap.isotopes.length; i++) {
			euclideanDistanceIntensities+=intensities[i]*intensities[i]; // add -1 into iDotP

			byte isotope=PrecursorScanMap.isotopes[i];
			if (isotope>=0) {
				// intensities[i] contains an extra -1 isotope
				isotopeDotProduct+=intensities[i]*predictedIsotopeDistribution[isotope];
				euclideanDistancePredicted+=predictedIsotopeDistribution[isotope]*predictedIsotopeDistribution[isotope];
			}
		}
		if (euclideanDistanceIntensities>0.0f&&euclideanDistancePredicted>0.0f) {
			euclideanDistanceIntensities=(float)Math.sqrt(euclideanDistanceIntensities);
			euclideanDistancePredicted=(float)Math.sqrt(euclideanDistancePredicted);
			isotopeDotProduct=isotopeDotProduct/(euclideanDistanceIntensities*euclideanDistancePredicted);
		} else {
			isotopeDotProduct=0.0f;
		}

		return new float[] {isotopeDotProduct, averagePPM};
	}
}
