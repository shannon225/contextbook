package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;

import org.jfree.chart.ChartPanel;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.CosineGaussian;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class EncyclopediaTwoScoringTask extends AbstractLibraryScoringTask {
	private final float dutyCycle;
	private final Range precursorIsolationRange;
	
	public EncyclopediaTwoScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue,
			SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.dutyCycle=dutyCycle;
		this.precursorIsolationRange=new Range(precursorIsolationRange.getStart(), precursorIsolationRange.getStop());
	}
	
	private static final int peaksKept=5;
	private static final int peaksConsidered=7;

	@Override
	protected Nothing process() {
		MassTolerance libraryTolerance=parameters.getLibraryFragmentTolerance();
		EncyclopediaScorer eScorer=(EncyclopediaScorer)scorer;
		int movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		for (LibraryEntry entry : super.entries) {
			if (parameters.getTopNTargetsUsed()>0) {
				entry=entry.trimToNPeaks(parameters.getTopNTargetsUsed(), parameters.getAAConstants());
			}
			
			AuxillaryPSMScorer auxScorer=eScorer.getAuxScorer().getEntryOptimizedScorer(entry);
			FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
			FragmentIon[] ions=model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true);
			Optional<FragmentIon[]> modificationSpecificIons;
			if (parameters.isVerifyModificationIons()) {
				modificationSpecificIons=model.getModificationSpecificIonObjects(precursorIsolationRange, parameters.getFragType(), entry.getPrecursorCharge(), true);
			} else {
				modificationSpecificIons=Optional.empty();
			}
			
			ions=FragmentIon.getUniqueFragments(ions, parameters.getFragmentTolerance()); // ensure that all ions are unique within tolerance

			double[] predictedMasses=entry.getMassArray();
			float[] predictedIntensities=entry.getIntensityArray();
			float[] correlation=entry.getCorrelationArray();
			ArrayList<FragmentIon> predictedTargets=new ArrayList<FragmentIon>();
			TFloatArrayList predictedTargetIntensities=new TFloatArrayList();
			TFloatArrayList maxCorrelations=new TFloatArrayList();
			for (FragmentIon target : ions) {
				int[] predictedIndicies=libraryTolerance.getIndicies(predictedMasses, target.getMass());
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
					maxCorrelations.add(maxCorrelation);
					predictedTargetIntensities.add(predictedIntensity);
					predictedTargets.add(target);
				}
			}
			if (predictedTargets.size()==0) continue;
			
			// limit to just measured in library
			ions=predictedTargets.toArray(new FragmentIon[0]);
			predictedIntensities=predictedTargetIntensities.toArray();
			correlation=maxCorrelations.toArray();
			double[] targetMasses=new double[ions.length];
			for (int i = 0; i < targetMasses.length; i++) {
				targetMasses[i]=ions[i].getMass();
			}
			
			AbstractScoringResult result=new PeptideScoringResult(entry);
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
			
			float[] primary=new float[super.stripes.size()];
			float[] secondary=new float[super.stripes.size()];
			double[][] allMasses=new double[super.stripes.size()][];
			float[][] allRawIntensities=new float[super.stripes.size()][];
			float[] retentionTimes=new float[super.stripes.size()];
			for (int i=0; i<super.stripes.size(); i++) {
				FragmentScan stripe=super.stripes.get(i);
				
				Triplet<float[], double[], float[]> results=score(stripe, ions, predictedIntensities, correlation);
				float[] scores=results.x;
				double[] masses=results.y;
				float[] intensities=results.z;
				
				primary[i]=scores[0];//*0.075f+scores[1]*0.43f;
				secondary[i]=scores[1];
				allMasses[i]=masses;
				allRawIntensities[i]=General.multiply(intensities, intensities); // undo the sqrt
				retentionTimes[i]=super.stripes.get(i).getScanStartTime();
				
				if (modificationSpecificIons.isPresent()) {
					// if modified signal represents less than 25% of the score then don't trust it
					float scoreFromModIons=eScorer.score(entry, stripe, modificationSpecificIons.get());
					if (scoreFromModIons/primary[i]<0.25f) {
						primary[i]=0.0f;
					}
				}
			}

			TFloatFloatHashMap map=new TFloatFloatHashMap();
			ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
			for (int i=0; i<primary.length; i++) {
				goodStripes.add(new ScoredIndex(primary[i], i));
				map.put(i, primary[i]);
			}
			Collections.sort(goodStripes);

			EValueCalculator calculator=new EValueCalculator(map, 0f, 0.5f);

			float bestScore=-Float.MAX_VALUE;
			float[] bestAuxScores=null;
			FragmentScan bestStripe=null;
			
			TIntHashSet takenScans=new TIntHashSet();
			int identifiedPeaks=0;
			int consideredPeaks=0;
			for (int i=goodStripes.size()-1; i>=0; i--) {
				float score=goodStripes.get(i).x;
				int index=goodStripes.get(i).y;
				if (takenScans.contains(index)) {
					continue;
					
				} else {	
					FragmentScan stripe=super.stripes.get(index);
					int lowerWindow=Math.max(0, index-2*movingAverageLength);
					int upperWindow=Math.min(super.stripes.size()-1, index+2*movingAverageLength);
					int range=Math.max(1, upperWindow-lowerWindow);
					float[] localRetentionTimes=new float[range];
					float[][] localRawIntensities=new float[range][]; // raw intensities, not sqrt
					System.arraycopy(retentionTimes, lowerWindow, localRetentionTimes, 0, range);
					System.arraycopy(allRawIntensities, lowerWindow, localRawIntensities, 0, range);
					
					float[][] chromatograms=General.transposeMatrix(localRawIntensities);
					ArrayList<float[]> chromatogramList=new ArrayList<float[]>();
					for (int j = 0; j < chromatograms.length; j++) {
						if (General.sum(chromatograms[j])>0.0f) {
							chromatograms[j]=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatograms[j]);
							if (parameters.isSubtractBackground()) {
								chromatograms[j]=backgroundSubtractMovingMedian(chromatograms[j], movingAverageLength*10);
							}
							chromatogramList.add(chromatograms[j]);
						}
					}
					
					TransitionRefinementData data=TransitionRefiner.identifyTransitions(entry.getPeptideModSeq(), entry.getPrecursorCharge(), stripe.getScanStartTime(), 
							ions, chromatogramList, localRetentionTimes, false, parameters);
					float[] correlations=data.getCorrelationArray();
					float[] integrations=data.getIntegrationArray();
					float[] median=data.getMedianChromatogram();
					Range rtRange=data.getRange();

//					System.out.println(stripe.getScanStartTime()/60+" --> "+rtRange.getStart()/60f+"\t"+rtRange.getStop()/60f);
//					HashMap<String, ChartPanel> panels=TransitionRefiner.getChartPanels(data);
//					Charter.launchCharts("TITLE", panels);
//					Charter.launchChart("RT (Min)", "Score", true, new XYTrace(General.divide(retentionTimes, 60f), primary, GraphType.line, "Primary"));					
//					Charter.launchChart("RT (Min)", "Intensity", true, new XYTrace(General.divide(localRetentionTimes, 60f), median, GraphType.line, "Median"));
//					try {
//						Thread.sleep(1000000);
//					} catch (Exception e) {}
					
					// 4 stdevs cover the range of 95% of the peak, assume peaks can't be smaller than 67% of the expected peak width
					float stdev = Math.max(parameters.getExpectedPeakWidth()/6.0f, rtRange.getRange()/4.0f);
					CosineGaussian gaussian=new CosineGaussian(stripe.getScanStartTime(), stdev, 1.0f);
					float[] idealGaussian=new float[median.length];
					float[] precursor=new float[median.length];
					float[] precursorPlusOne=new float[median.length];
					float sumMedianInRange=0;
					for (int j = 0; j < idealGaussian.length; j++) {
						idealGaussian[j]=(float)gaussian.getPDF(localRetentionTimes[j]);
						if (rtRange.contains(localRetentionTimes[j])) {
							sumMedianInRange+=median[j];
							
							// zero outside peak
							Peak[] peaks=precursors.getIsotopePacket(entry.getPrecursorMZ(), localRetentionTimes[j], entry.getPrecursorCharge(), parameters.getPrecursorTolerance());
							precursor[j]=peaks[PrecursorScanMap.monoisotopicIndex].intensity;
							precursorPlusOne[j]=peaks[PrecursorScanMap.plusOneIndex].intensity;
						}
					}
					
					//Spectrum acquiredPeak=getSpectrum(allMasses[index], integrations, stripe.getScanStartTime(), stripe.getPrecursorMZ());
					
					float[] auxScoreArray=auxScorer.score(entry, stripe, predictedIsotopeDistribution, precursors);
					
					// ADDITIONAL AUX SCORES (make sure to add titles)
					
					// main score evalue
					float evalue=calculator.getNegLnEValue(score);
					if (Float.isNaN(evalue)) {
						evalue=-1.0f;
					}
					
					// correlation scores
					float[] sortedCorrelations=correlations.clone();
					Arrays.sort(sortedCorrelations);
					General.reverse(sortedCorrelations);
					
					int numPeaksWithGoodCorrelation=0;
					int numPeaksWithGreatCorrelation=0;
					for (int j = 0; j < sortedCorrelations.length; j++) {
						if (sortedCorrelations[j]>=TransitionRefiner.identificationCorrelationThreshold) {
							numPeaksWithGoodCorrelation++;
						}

						if (sortedCorrelations[j]>=TransitionRefiner.quantitativeCorrelationThreshold) {
							numPeaksWithGreatCorrelation++;
						}
					}
					

					int isIntegratedSignal=0;
					int isIntegratedPrecursor=0;
					float correlationToGaussian=0.0f;
					float correlationToPrecursor=0.0f;
					//float correlationToPlusOne=0.0f;
					if (sumMedianInRange>0) {
						if (numPeaksWithGreatCorrelation>0) {
							isIntegratedSignal=1;
						}
						if (General.sum(precursor)>0) {
							isIntegratedPrecursor=1;
						}
						correlationToGaussian=Correlation.getPearsons(idealGaussian, median);
						correlationToPrecursor=Correlation.getPearsons(precursor, median);
						//correlationToPlusOne=Correlation.getPearsons(precursorPlusOne, median);
					}
					
					auxScoreArray=General.concatenate(auxScoreArray, new float[] {evalue, correlationToGaussian, 
							correlationToPrecursor, isIntegratedSignal, isIntegratedPrecursor,
							numPeaksWithGoodCorrelation});
					
					// block out a 40 scan window
					for (int j=lowerWindow; j<=upperWindow; j++) {
						takenScans.add(j);
					}
					if (bestStripe==null) {
						bestStripe=stripe;
						bestAuxScores=auxScoreArray;
						bestScore=score;
					}

					if (consideredPeaks>peaksConsidered) {
						break;
					}
					consideredPeaks++;
					
					// don't add if we don't have enough quantitative ions (add below if it's still the best and there are no matches)
					if (numPeaksWithGoodCorrelation<parameters.getMinNumOfQuantitativePeaks()) continue;
					
					float deltaPrecursorMass=auxScorer.getParentDeltaMassIndex()>=0?auxScoreArray[auxScorer.getParentDeltaMassIndex()]:0.0f;
					float deltaFragmentMass=auxScorer.getFragmentDeltaMassIndex()>=0?auxScoreArray[auxScorer.getFragmentDeltaMassIndex()]:0.0f;
					result.addStripe(score, auxScoreArray, deltaPrecursorMass, deltaFragmentMass, stripe);
					if (identifiedPeaks>peaksKept) {
						// keep N+1 peaks
						break;
					}
					identifiedPeaks++;
				}
			}
			
			if (identifiedPeaks==0) {
				// add the best data if we can't find any valid peaks
				float deltaPrecursorMass=auxScorer.getParentDeltaMassIndex()>=0?bestAuxScores[auxScorer.getParentDeltaMassIndex()]:0.0f;
				float deltaFragmentMass=auxScorer.getFragmentDeltaMassIndex()>=0?bestAuxScores[auxScorer.getFragmentDeltaMassIndex()]:0.0f;
				result.addStripe(bestScore, bestAuxScores, deltaPrecursorMass, deltaFragmentMass, bestStripe);
			}
			
			resultsQueue.add(result);
		}
		return Nothing.NOTHING;
	}

//	public TransitionRefinementData pickPeak(LibraryEntry unitEntry) {
//		// invert each scan into fragment chromatograms (RTs by fragment)
//		ArrayList<PeakScores> bestKeptPeaks=new ArrayList<PeakScores>();
//		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
//		ArrayList<FragmentIon> fragmentMasses=new ArrayList<FragmentIon>();
//		for (int i=0; i<bestScores.length; i++) {
//			if (bestScores[i]!=null&&bestScores[i].getScore()>0) {
//				float[] chromatogram=traces[i].toArray();
//				chromatogram=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatogram);
//				chromatograms.add(chromatogram);
//				bestKeptPeaks.add(bestScores[i]);
//				fragmentMasses.add(bestScores[i].getTarget());
//			}
//		}
//
//		// identify transitions
//		TransitionRefinementData data=TransitionRefiner.identifyTransitions(unitEntry.getPeptideModSeq(), unitEntry.getPrecursorCharge(), unitEntry.getScanStartTime(), fragmentMasses.toArray(new FragmentIon[fragmentMasses.size()]), chromatograms, retentionTimes.toArray(), wasInferred, params);
//		float[] correlations=data.getCorrelationArray();
//		float[] integrations=data.getIntegrationArray();
//		Range rtRange=data.getRange();
//	}

	public Triplet<float[], double[], float[]> score(Spectrum spectrum, FragmentIon[] ions, float[] predictedTargetIntensities, float[] correlations) {
		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		TDoubleArrayList actualTargetMasses=new TDoubleArrayList();
		TFloatArrayList actualTargetIntensities=new TFloatArrayList();
		for (int i = 0; i < ions.length; i++) {
			FragmentIon target=ions[i];
		
			int[] indicies=acquiredTolerance.getIndicies(acquiredMasses, target.getMass());
			float intensity=0.0f;
			float bestPeakIntensity=0.0f;
			double bestPeakMass=0.0;
			
			for (int j=0; j<indicies.length; j++) {
				intensity+=acquiredIntensities[indicies[j]];
				
				if (acquiredIntensities[indicies[j]]>bestPeakIntensity) {
					bestPeakIntensity=acquiredIntensities[indicies[j]];
					bestPeakMass=acquiredMasses[indicies[j]];
				}
			}
			actualTargetIntensities.add(intensity);
			actualTargetMasses.add(bestPeakMass);
		}

		double[] actualTargetMassesRaw=actualTargetMasses.toArray();
		float[] actualTargetIntensitiesRaw=actualTargetIntensities.toArray();
		
		float[] predictedTargetIntensitiesArray=General.normalizeToL2(predictedTargetIntensities);
		float[] actualTargetIntensitiesArray=General.normalizeToL2(actualTargetIntensitiesRaw);
		
		float dotProduct=General.sum(General.multiply(General.multiply(predictedTargetIntensitiesArray, actualTargetIntensitiesArray), correlations));

		if (Float.isNaN(dotProduct)||dotProduct<0.0f) dotProduct=0.0f;
		
		float sumOfSquaredErrors=0.0f; // normalized to sum of targeted intensities

		int numberOfMatchingPeaks=0;
		for (int i=0; i<predictedTargetIntensitiesArray.length; i++) {
			if (predictedTargetIntensitiesArray[i]>0.0) {
				float delta=predictedTargetIntensitiesArray[i]-actualTargetIntensitiesArray[i];
				float deltaSquared=delta*delta;
				sumOfSquaredErrors+=deltaSquared;
				if (actualTargetIntensitiesArray[i]>0.0) {
					numberOfMatchingPeaks++;
				}
			}
		}
		
		float xTandem;
		if (numberOfMatchingPeaks==0||dotProduct<=0) {
			xTandem=0.0f;
		} else {
			xTandem=((float)Log.protectedLog10(dotProduct))+Log.logFactorial(numberOfMatchingPeaks); // really log10(X!Tandem score)
		}
		
		float scribe;
		if (sumOfSquaredErrors<=0.0f) {
			scribe=0.0f;
		} else {
			scribe=Log.protectedLn(1.0f/sumOfSquaredErrors);
		}
		
		return new Triplet<float[], double[], float[]>(new float[] {xTandem, scribe}, 
				actualTargetMassesRaw, actualTargetIntensitiesRaw);
	}


	static Spectrum getSpectrum(final double[] masses, final float[] intensities, final float scanStartTime, final double mz) {
		return new Spectrum() {
			@Override
			public float getTIC() {
				return General.sum(intensities);
			}

			@Override
			public String getSpectrumName() {
				return Double.toString(mz)+","+Float.toString(scanStartTime);
			}

			@Override
			public float getScanStartTime() {
				return scanStartTime;
			}

			@Override
			public double getPrecursorMZ() {
				return mz;
			}

			@Override
			public double[] getMassArray() {
				return masses;
			}

			@Override
			public float[] getIntensityArray() {
				return intensities;
			}
		};
	}
}
