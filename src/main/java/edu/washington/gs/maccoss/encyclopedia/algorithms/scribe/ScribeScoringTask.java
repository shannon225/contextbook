package edu.washington.gs.maccoss.encyclopedia.algorithms.scribe;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.Scribe;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.SimilarPeptideBinner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class ScribeScoringTask extends AbstractLibraryScoringTask {
	// FIXME shouldn't be static (long term)
	private static final HashMap<String, SparseXCorrSpectrum> libraryEntryMap=new HashMap<String, SparseXCorrSpectrum>();
	private static final HashMap<String, float[]> isotopeDistributions=new HashMap<String, float[]>();
	private final ScribeScorer scorerFunction;
	private final Range fragmentRange=new Range(200, Float.MAX_VALUE); // remove peaks below 200 
	private static final int minimumNumberOfPeaks = 10;
	
	private static final int peaksKept=1;
	
	public ScribeScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue, SearchParameters parameters) {
		super(scorer, filterEntriesByScore(entries), stripes, precursors, resultsQueue, parameters);
		scorerFunction=(ScribeScorer)scorer;
	}
	
	public static ArrayList<LibraryEntry> filterEntriesByScore(ArrayList<LibraryEntry> entries) {
		HashMap<String, LibraryEntry> map=new HashMap<>();
		for (LibraryEntry entry : entries) {
			LibraryEntry prev=map.get(entry.getPeptideModSeq());
			if (prev==null||prev.getTIC()<entry.getTIC()) {
				// TODO Smart to prefer higher TIC? at least this keeps a consistent sort order
				map.put(entry.getPeptideModSeq(), entry);
			}
		}
		return new ArrayList<LibraryEntry>(map.values());
	}

	@Override
	protected Nothing process() {
		for (FragmentScan msms : super.stripes) {
			msms=msms.trimMasses(fragmentRange);
			XCorrCalculatorSpectrum xcorrMSMS=new XCorrCalculatorSpectrum(msms, parameters);
			if (xcorrMSMS.getNumPeaks()<minimumNumberOfPeaks) {
				//System.out.println("cut\t"+1);
				continue;
			}
			
			ArrayList<ScoredIndex> goodHits=new ArrayList<ScoredIndex>();
			TFloatFloatHashMap map=new TFloatFloatHashMap(); // x=index, y=score
			float maxXCorr=0.0f;
			String maxSequence=null;
			float secondMaxXCorr=0.0f;
			float[] xcorrs=new float[super.entries.size()];
			for (int i=0; i<super.entries.size(); i++) {
				LibraryEntry entry=super.entries.get(i);
				boolean match=parameters.getPrecursorTolerance().equals(entry.getPrecursorMZ(), msms.getPrecursorMZ());
				for (int j = 0; j < Scribe.NUMBER_OF_ISOTOPES_ABOVE_MONOISOTOPIC; j++) {
					double target = entry.getPrecursorMZ()+(j+1)*MassConstants.neutronMass/entry.getPrecursorCharge();
					match=match||parameters.getPrecursorTolerance().equals(target, msms.getPrecursorMZ());	
				}
				if (match) {
					SparseXCorrSpectrum xcorrEntry=getXCorrEntry(entry);
					
					float score=xcorrMSMS.score(xcorrEntry);
					xcorrs[i]=score;
					float[] otherScores=score(entry, msms);
					
					if (otherScores[0]>0) {
						float composite=otherScores[1];
						goodHits.add(new ScoredIndex(composite, i));
						map.put(i, composite);
					}
					
					if (maxSequence==null) {
						maxSequence=entry.getPeptideSeq();
						maxXCorr=score;
					} else if (score>maxXCorr) {
						if (!SimilarPeptideBinner.areSimilarEnough(maxSequence, entry.getPeptideSeq())) {
							secondMaxXCorr=maxXCorr;
							maxSequence=entry.getPeptideSeq();
							maxXCorr=score;
						}
					} else if (score>secondMaxXCorr) {
						if (!SimilarPeptideBinner.areSimilarEnough(maxSequence, entry.getPeptideSeq())) {
							secondMaxXCorr=score;
						}
					}
				}
			}
			if (map.size()==0) {
				//System.out.println("cut\t"+2);
				continue;
			}
			
			EValueCalculator calculator=new EValueCalculator(map, 0.1f, 0.1f);
			
			Collections.sort(goodHits);
			int identifiedPeaks=0;
			for (int i=goodHits.size()-1; i>=0; i--) {
				float score=goodHits.get(i).x;
				int index=goodHits.get(i).y;
				float evalue=calculator.getNegLnEValue(score);
				float xcorr=xcorrs[index];
				
				float deltaCn=(xcorr==0.0f||secondMaxXCorr==0.0f)?0.0f:(xcorr-secondMaxXCorr)/maxXCorr;

				LibraryEntry entry=super.entries.get(index);
					
				float[] predictedIsotopeDistribution=getIsotopeDistribution(entry);
				float[] auxScoreArray=scorerFunction.auxScore(entry, msms, predictedIsotopeDistribution, precursors);
				if (auxScoreArray[0]<=0) {
					// dot product is 0, means no matching b/y peaks in the top scoring model
					//System.out.println("cut\t"+3);
					//continue;
				}

				AbstractScoringResult result=new PeptideScoringResult(entry);
				result.addStripe(score, General.concatenate(auxScoreArray, evalue, map.size(), deltaCn), msms);
				
				if (identifiedPeaks>peaksKept) {
					// keep N+1 peaks
					break;
				}
				identifiedPeaks++;
				resultsQueue.add(result);
			}
		}
		return Nothing.NOTHING;
	}
	public SparseXCorrSpectrum getXCorrEntry(LibraryEntry entry) {
		SparseXCorrSpectrum xcorr=libraryEntryMap.get(entry.getPeptideModSeq());
		if (xcorr==null) {
			xcorr=SparseXCorrCalculator.normalize(entry, new Range((float)entry.getPrecursorMZ()-10f, (float)entry.getPrecursorMZ()+10f), false, parameters);
			libraryEntryMap.put(entry.getPeptideModSeq(), xcorr);
		}
		return xcorr;
		
	}

	public float[] getIsotopeDistribution(LibraryEntry entry) {
		float[] predictedIsotopeDistribution=isotopeDistributions.get(entry.getPeptideModSeq());
		if (predictedIsotopeDistribution==null) {
			predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());
			isotopeDistributions.put(entry.getPeptideModSeq(), predictedIsotopeDistribution);
		}
		return predictedIsotopeDistribution;
	}


	public float[] score(LibraryEntry entry, Spectrum spectrum) {
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
				for (int j=0; j<indicies.length; j++) {
					intensity+=acquiredIntensities[indicies[j]];
					
					if (acquiredIntensities[indicies[j]]>bestPeakIntensity) {
						bestPeakIntensity=acquiredIntensities[indicies[j]];
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
			}
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
		
		if (sumOfSquaredErrors<1e-5f) {
			sumOfSquaredErrors=1e-5f;
		}

		float xTandem;
		if (numberOfMatchingPeaks==0) {
			xTandem=0.0f;
		} else {
			xTandem=((float)Log.protectedLog10(dotProduct))+Log.logFactorial(numberOfMatchingPeaks); // really log10(X!Tandem score)
		}
		
		return new float[] {xTandem, Log.protectedLn(1.0f/sumOfSquaredErrors)};
	}
}
