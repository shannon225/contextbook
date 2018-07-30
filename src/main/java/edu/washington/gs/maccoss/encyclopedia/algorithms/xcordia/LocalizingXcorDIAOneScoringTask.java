package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.MedianChromatogramData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.VariantFastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.IndexedObject;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class LocalizingXcorDIAOneScoringTask extends AbstractLibraryScoringTask {
	private static final int MINIMUM_NUMBER_OF_PEAKS = 3;
	private static final float MINIMUM_LOCALIZATION_SCORE_TO_CONTINUE = 2.0f;
	
	private final BackgroundFrequencyInterface background;
	private final BlockingQueue<ModificationLocalizationData> localizationQueue;
	private final int movingAverageLength;
	
	public LocalizingXcorDIAOneScoringTask(PSMScorer scorer, BackgroundFrequencyInterface background, ArrayList<LibraryEntry> entries, 
			ArrayList<Stripe> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue,
			BlockingQueue<ModificationLocalizationData> localizationQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.background=background;
		this.localizationQueue=localizationQueue;
		this.movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		
		assert(scorer instanceof XCorDIAOneScorer);
	}

	@Override
	protected Nothing process() {
		// separate targets from decoys and process in batches
		ArrayList<XCorrLibraryEntry> targetBatch=new ArrayList<>();
		ArrayList<XCorrLibraryEntry> decoyBatch=new ArrayList<>();
		for (LibraryEntry entry : super.entries) {
			XCorrLibraryEntry xcordiaEntry = getXCorrEntry(entry, parameters);
			if (entry.isDecoy()) {
				decoyBatch.add(xcordiaEntry);
			} else {
				targetBatch.add(xcordiaEntry);
			}
		}
		processPeptide(targetBatch);
		processPeptide(decoyBatch);

		return Nothing.NOTHING;
	}
	
	private void processPeptide(ArrayList<XCorrLibraryEntry> seedEntries) {
		// score N peptides individually
		ArrayList<ScoredObject<IndexedObject<PeptidePrecursor>>> bestScoresByEntry=new ArrayList<>();
		HashMap<PeptidePrecursor, float[]> scoresByEntry=new HashMap<>();
		HashMap<PeptidePrecursor, float[]> isotopesByEntry=new HashMap<>();
		HashMap<PeptidePrecursor, FragmentationModel> modelsByEntry=new HashMap<>();
		for (XCorrLibraryEntry xcordiaEntry : seedEntries) {
			modelsByEntry.put(xcordiaEntry, PeptideUtils.getPeptideModel(xcordiaEntry.getPeptideModSeq(), parameters.getAAConstants()));
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(xcordiaEntry.getPeptideModSeq(), parameters.getAAConstants());
			isotopesByEntry.put(xcordiaEntry, predictedIsotopeDistribution);
			
			// perform initial scoring
			float[] primary=new float[super.stripes.size()];
			for (int i=0; i<super.stripes.size(); i++) {
				Stripe stripe=super.stripes.get(i);
				XCorrStripe xcordiaStripe;
				if (stripe instanceof XCorrStripe) {
					xcordiaStripe=(XCorrStripe)stripe;
				} else {
					xcordiaStripe=new XCorrStripe(stripe, parameters);
				}
				primary[i]=scorer.score(xcordiaEntry, xcordiaStripe, predictedIsotopeDistribution, precursors);
			}
			
			float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);
			
			// determine the N best peaks for this peptide
			LinkedList<ScoredIndex> keptIndicies=new LinkedList<>();
			for (int i=0; i<averagePrimary.length; i++) {
				if (keptIndicies.size()==0) {
					keptIndicies.add(new ScoredIndex(averagePrimary[i], i));
				} else if (averagePrimary[i]>keptIndicies.get(keptIndicies.size()-1).x) {
					int upperIndexRange=i+movingAverageLength;
					int lowerIndexRange=i-movingAverageLength;
					int count=0;
					for (ScoredIndex scoredIndex : keptIndicies) {
						if (averagePrimary[i]>scoredIndex.x&&(scoredIndex.y<lowerIndexRange||scoredIndex.y>upperIndexRange)) {
							break;
						}
						count++;
					}
					keptIndicies.add(count, new ScoredIndex(averagePrimary[i], i));
					if (keptIndicies.size()>seedEntries.size()) {
						keptIndicies.removeLast();
					}
				}
			}
			
			for (ScoredIndex scoredIndex : keptIndicies) {
				bestScoresByEntry.add(new ScoredObject<IndexedObject<PeptidePrecursor>>(scoredIndex.x, new IndexedObject<PeptidePrecursor>(scoredIndex.y, xcordiaEntry)));
			}
			scoresByEntry.put(xcordiaEntry, averagePrimary);
		}
		Collections.sort(bestScoresByEntry);
		
		// if there's only one peptide then no need to compare against anything
		if (seedEntries.size()==1) {
			XCorrLibraryEntry xcordiaEntry=(XCorrLibraryEntry)seedEntries.get(0);
			PeptideScoringResult result=new PeptideScoringResult(xcordiaEntry);
			
			float[] xcorrArray = scoresByEntry.get(xcordiaEntry);
			TFloatFloatHashMap map=new TFloatFloatHashMap();
			float maxXCorr=-Float.MAX_VALUE;
			int maxXCorrIndex=0;
			for (int i=0; i<xcorrArray.length; i++) {
				if (xcorrArray[i]>maxXCorr) {
					maxXCorr=xcorrArray[i];
					maxXCorrIndex=i;
				}
				map.put(i, xcorrArray[i]);
			}
			EValueCalculator evalueCalculator=new EValueCalculator(map);
			float evalue=evalueCalculator.getNegLog10EValue(maxXCorr);
			Stripe stripe = stripes.get(maxXCorrIndex);
			float[] auxScoreArray=scorer.auxScore(xcordiaEntry, stripe, isotopesByEntry.get(xcordiaEntry), precursors);

			FragmentationModel topModel=modelsByEntry.get(xcordiaEntry);
			FragmentIon[] allIonObjectsWithMissing = topModel.getPrimaryIonObjects(parameters.getFragType(), xcordiaEntry.getPrecursorCharge(), true);
			double[] allIonsWithMissing = FragmentIon.getMasses(allIonObjectsWithMissing);
			ChromatogramMap chromMap=new ChromatogramMap(allIonsWithMissing, maxXCorrIndex);
			ArrayList<float[]> allChromatogramsWithMissing=chromMap.getChromatograms(allIonsWithMissing);
			
			ArrayList<FragmentIon> allIonObjects=new ArrayList<>();
			ArrayList<float[]> allChromatograms=new ArrayList<>();
			for (int i = 0; i < allIonsWithMissing.length; i++) {
				float[] chromatogram = allChromatogramsWithMissing.get(i);
				if (testChromatogramForValidity(chromatogram)) {
					allIonObjects.add(allIonObjectsWithMissing[i]);
					allChromatograms.add(chromatogram);
				}
			}
			MedianChromatogramData allIonsData=TransitionRefiner.extractMedianChromatogram(stripe.getScanStartTime(), allChromatograms, chromMap.getRetentionTimes(), Optional.empty(), true);
			float medianMean=General.mean(allIonsData.getMedianChromatogram(), allIonsData.getIndices().getStart(), allIonsData.getIndices().getStop());
			int numberOfPeaks=0;
			float sumCorrelation=0.0f;
			ArrayList<FragmentIon> foundLocalizingIons=new ArrayList<>();
			for (int i=0; i<allIonsData.getNormalizedChromatograms().size(); i++) {
				float[] normalizedChromatogram=allIonsData.getNormalizedChromatograms().get(i);
				float correlation=TransitionRefiner.calculateCorrelation(medianMean, allIonsData.getIndices(), allIonsData.getMedianChromatogram(), normalizedChromatogram);
				if (correlation>TransitionRefiner.identificationCorrelationThreshold) {
					numberOfPeaks++;
					foundLocalizingIons.add(allIonObjects.get(i));
				}
				sumCorrelation+=correlation*correlation;
			}
			
			result.addStripe(maxXCorr, General.concatenate(auxScoreArray, maxXCorr, evalue, seedEntries.size()>1?1:0, numberOfPeaks, sumCorrelation), stripe);
			resultsQueue.add(result);

			String annotation;
			PeptideModification localizingModification;
			if (xcordiaEntry.getPeptide() instanceof VariantFastaPeptideEntry) {
				localizingModification=PeptideModification.polymorphism; 
				annotation=((VariantFastaPeptideEntry)xcordiaEntry.getPeptide()).getVariant().toString();
			} else if (parameters.getLocalizingModification().isPresent()) {
				localizingModification=parameters.getLocalizingModification().get();
				annotation=localizingModification.getName();
			} else {
				localizingModification=PeptideModification.polymorphism; 
				annotation="canonical";
			}
			AmbiguousPeptideModSeq ambiguousPeptideModSeq=AmbiguousPeptideModSeq.getUnambigous(xcordiaEntry.getPeptideModSeq(), localizingModification, parameters.getAAConstants(), annotation);
			boolean isSiteSpecific = true;
			ModificationLocalizationData modData=new ModificationLocalizationData(ambiguousPeptideModSeq, stripe.getScanStartTime(), 1000.0f, numberOfPeaks, 0, isSiteSpecific, isSiteSpecific, !isSiteSpecific, foundLocalizingIons.toArray(new FragmentIon[foundLocalizingIons.size()]), 0.0f, 0.0f);
			try {
				localizationQueue.put(modData);
			} catch (InterruptedException e) {
				Logger.logException(e);
			}
			return; // return out to skip localization processing
		}
		
		// determine N total peaks worth scoring (at least one for each peptide)
		TIntArrayList pickedPeakIndicies=new TIntArrayList();
		HashSet<PeptidePrecursor> previouslyPicked=new HashSet<>();
		for (int i=bestScoresByEntry.size()-1; i>=0; i--) {
			ScoredObject<IndexedObject<PeptidePrecursor>> scoredPeptide=bestScoresByEntry.get(i);
			
			if (!previouslyPicked.contains(scoredPeptide.y.y)) {
				previouslyPicked.add(scoredPeptide.y.y);
				pickedPeakIndicies.add(scoredPeptide.y.x);
			}
			
			if (previouslyPicked.size()==seedEntries.size()) break;
		}

		// localize each of N peaks
		float bestXCorrSoFar=-Float.MAX_VALUE;
		PeptidePrecursor bestPeptideSoFar=null;
		ScoreData bestScoreSoFar=null;
		HashMap<PeptidePrecursor, ScoreData> scoredPeptides=new HashMap<>();
		for (int index : pickedPeakIndicies.toArray()) {
			// find the top and second best scoring peptides
			float topScore=-Float.MAX_VALUE;
			float secondTopScore=-Float.MAX_VALUE;
			PeptidePrecursor topScoring=null;
			PeptidePrecursor secondTopScoring=null;

			for (XCorrLibraryEntry xcordiaEntry : seedEntries) {
				float score=scoresByEntry.get(xcordiaEntry)[index];
				if (score>topScore) {
					secondTopScore=topScore;
					secondTopScoring=topScoring;
					topScore=score;
					topScoring=xcordiaEntry;
				} else if (score>secondTopScore) {
					secondTopScore=score;
					secondTopScoring=xcordiaEntry;
				}
			}
			
			if (scoredPeptides.containsKey(topScoring)) continue;

			// get fragment ion map containing all ions
			FragmentationModel topModel=modelsByEntry.get(topScoring);
			FragmentationModel secondTopModel=modelsByEntry.get(secondTopScoring);
			FragmentIon[] allIonObjectsWithMissing = topModel.getPrimaryIonObjects(parameters.getFragType(), topScoring.getPrecursorCharge(), true);
			double[] allIonsWithMissing = FragmentIon.getMasses(allIonObjectsWithMissing);
			ChromatogramMap chromMap=new ChromatogramMap(allIonsWithMissing, index);
			
			// variant-specific median trace
			FragmentIon[] uniqueIonsWithMissing=getUniqueFragmentIons(topModel, secondTopModel, topScoring.getPrecursorCharge(), parameters.getFragType());
			double[] uniqueIonMassesWithMissing = FragmentIon.getMasses(uniqueIonsWithMissing);
			ArrayList<float[]> variantSpecificChromatogramsWithMissing=chromMap.getChromatograms(uniqueIonMassesWithMissing);

			ArrayList<FragmentIon> uniqueIons=new ArrayList<>();
			ArrayList<float[]> variantSpecificChromatograms=new ArrayList<>();
			for (int i = 0; i < uniqueIonsWithMissing.length; i++) {
				float[] chromatogram = variantSpecificChromatogramsWithMissing.get(i);
				if (testChromatogramForValidity(chromatogram)) {
					uniqueIons.add(uniqueIonsWithMissing[i]);
					variantSpecificChromatograms.add(chromatogram);
				}
			}
			double[] uniqueIonMasses=FragmentIon.getMasses(uniqueIons.toArray(new FragmentIon[uniqueIons.size()]));
			
			MedianChromatogramData variantSpecificData=TransitionRefiner.extractMedianChromatogram(stripes.get(index).getScanStartTime(), variantSpecificChromatograms, chromMap.getRetentionTimes(), Optional.empty(), true);
			float medianMean=General.mean(variantSpecificData.getMedianChromatogram(), variantSpecificData.getIndices().getStart(), variantSpecificData.getIndices().getStop());

			// calculate variant-specific p-value
			float[] frequencies=background.getFrequencies(uniqueIonMasses, ((XCorrLibraryEntry)topScoring).getPrecursorMZ(), parameters.getFragmentTolerance());
			float logProb=0;
			ArrayList<FragmentIon> foundLocalizingIons=new ArrayList<>(); 
			if (medianMean>0) {
				for (int j = 0; j < uniqueIonMasses.length; j++) {
					float[] normalizedChromatogram=variantSpecificData.getNormalizedChromatograms().get(j);
					float correlation=TransitionRefiner.calculateCorrelation(medianMean, variantSpecificData.getIndices(), variantSpecificData.getMedianChromatogram(), normalizedChromatogram);
					if (correlation>TransitionRefiner.identificationCorrelationThreshold) {
						logProb+=-Log.protectedLog10(frequencies[j]);
						foundLocalizingIons.add(uniqueIons.get(j));
					}
				}
			}

			if (logProb>=MINIMUM_LOCALIZATION_SCORE_TO_CONTINUE||topScore>bestXCorrSoFar) {
				// check all ions versus median trace
				ArrayList<float[]> allChromatogramsWithMissing=chromMap.getChromatograms(allIonsWithMissing);
				
				ArrayList<float[]> allChromatograms=new ArrayList<>();
				for (int i = 0; i < allIonsWithMissing.length; i++) {
					float[] chromatogram = allChromatogramsWithMissing.get(i);
					if (testChromatogramForValidity(chromatogram)) {
						allChromatograms.add(chromatogram);
					}
				}
				Optional<float[]> siteSpecificPeakshape;
				if (logProb<MINIMUM_LOCALIZATION_SCORE_TO_CONTINUE) {
					// if we don't trust the localization, just fit the peak to all the ions
					siteSpecificPeakshape = Optional.empty();
				} else {
					siteSpecificPeakshape = Optional.of(variantSpecificData.getMedianChromatogram());
				}
				
				MedianChromatogramData allIonsData=TransitionRefiner.extractMedianChromatogram(stripes.get(index).getScanStartTime(), allChromatograms, chromMap.getRetentionTimes(), siteSpecificPeakshape, true);
				medianMean=General.mean(allIonsData.getMedianChromatogram(), allIonsData.getIndices().getStart(), allIonsData.getIndices().getStop());
				int numberOfPeaks=0;
				float sumCorrelation=0.0f;
				for (int i=0; i<allIonsData.getNormalizedChromatograms().size(); i++) {
					float[] normalizedChromatogram=allIonsData.getNormalizedChromatograms().get(i);
					float correlation=TransitionRefiner.calculateCorrelation(medianMean, allIonsData.getIndices(), allIonsData.getMedianChromatogram(), normalizedChromatogram);
					if (correlation>TransitionRefiner.identificationCorrelationThreshold) {
						numberOfPeaks++;
					}
					sumCorrelation+=correlation*correlation;
				}
				
				ScoreData score=new ScoreData(logProb, topScore, index, numberOfPeaks, sumCorrelation, foundLocalizingIons);
				// independently keep track of the best peptide, even if it doesn't pass localization thresholds
				if (topScore>bestXCorrSoFar) {
					bestPeptideSoFar=topScoring;
					bestScoreSoFar=score;
				}
				if (logProb>=MINIMUM_LOCALIZATION_SCORE_TO_CONTINUE&&numberOfPeaks>=MINIMUM_NUMBER_OF_PEAKS) {
					scoredPeptides.put(topScoring, score);
				}
			}
		}
		
		// make sure at least one peptide is scored, even if it doesn't pass the thresholds
		boolean needToAddUnlocalizedPeptide=scoredPeptides.size()==0&&bestPeptideSoFar!=null;
		if (needToAddUnlocalizedPeptide) { 
			scoredPeptides.put(bestPeptideSoFar,	bestScoreSoFar);
		}

		// process scored peptides
		for (Entry<PeptidePrecursor, ScoreData> entry : scoredPeptides.entrySet()) {
			XCorrLibraryEntry xcordiaEntry=(XCorrLibraryEntry) entry.getKey();
			ScoreData data=entry.getValue();
			PeptideScoringResult result=new PeptideScoringResult(xcordiaEntry);
			
			float[] xcorrArray = scoresByEntry.get(xcordiaEntry);
			TFloatFloatHashMap map=new TFloatFloatHashMap();
			for (int i=0; i<xcorrArray.length; i++) {
				map.put(i, xcorrArray[i]);
			}
			EValueCalculator evalueCalculator=new EValueCalculator(map);
			float evalue=evalueCalculator.getNegLog10EValue(data.xCorr);
			Stripe stripe = stripes.get(data.spectrumIndex);
			float[] auxScoreArray=scorer.auxScore(xcordiaEntry, stripe, isotopesByEntry.get(xcordiaEntry), precursors);
			
			result.addStripe(data.xCorr, General.concatenate(auxScoreArray, data.xCorr, evalue, seedEntries.size()>1?1:0, data.numberOfPeaks, data.sumCorrelation), stripe);
			resultsQueue.add(result);

			String annotation;
			PeptideModification localizingModification;
			if (xcordiaEntry.getPeptide() instanceof VariantFastaPeptideEntry) {
				localizingModification=PeptideModification.polymorphism; 
				annotation=((VariantFastaPeptideEntry)xcordiaEntry.getPeptide()).getVariant().toString();
			} else if (parameters.getLocalizingModification().isPresent()) {
				localizingModification=parameters.getLocalizingModification().get();
				annotation=localizingModification.getName();
			} else {
				localizingModification=PeptideModification.polymorphism; 
				annotation="canonical";
			}
			AmbiguousPeptideModSeq ambiguousPeptideModSeq=AmbiguousPeptideModSeq.getUnambigous(xcordiaEntry.getPeptideModSeq(), localizingModification, parameters.getAAConstants(), annotation);
			boolean isSiteSpecific = data.localizationScore>=MINIMUM_LOCALIZATION_SCORE_TO_CONTINUE&&data.numberOfPeaks>=MINIMUM_NUMBER_OF_PEAKS;
			ModificationLocalizationData modData=new ModificationLocalizationData(ambiguousPeptideModSeq, stripe.getScanStartTime(), data.localizationScore, data.numberOfPeaks, 0, isSiteSpecific, isSiteSpecific, !isSiteSpecific, data.foundLocalizingIons.toArray(new FragmentIon[data.foundLocalizingIons.size()]), 0.0f, 0.0f);
			try {
				localizationQueue.put(modData);
			} catch (InterruptedException e) {
				Logger.logException(e);
			}
		}
	}

	// require at least two consistent values
	private boolean testChromatogramForValidity(float[] chromatogram) {
		for (int i = 1; i < chromatogram.length; i++) {
			if (chromatogram[i-1]>0.0f&&chromatogram[i]>0.0f) {
				return true;
			}
		}
		return false;
	}
	
	public static FragmentIon[] getUniqueFragmentIons(FragmentationModel topModel, FragmentationModel secondTopModel, byte precursorCharge, FragmentationType fragType) {
		HashSet<FragmentIon> ions=new HashSet<FragmentIon>(Arrays.asList(topModel.getPrimaryIonObjects(fragType, precursorCharge, false, true)));
		ions.removeAll(Arrays.asList(secondTopModel.getPrimaryIonObjects(fragType, precursorCharge, false, true)));

		FragmentIon[] ionArray=ions.toArray(new FragmentIon[ions.size()]);
		Arrays.sort(ionArray);
		return ionArray;
	}
	
	private class ScoreData {
		final float localizationScore;
		final int spectrumIndex;
		final int numberOfPeaks;
		final float sumCorrelation;
		final float xCorr;
		final ArrayList<FragmentIon> foundLocalizingIons;

		public ScoreData(float localizationScore, float xCorr, int spectrumIndex, int numberOfPeaks, float sumCorrelation, ArrayList<FragmentIon> foundLocalizingIons) {
			this.localizationScore=localizationScore;
			this.xCorr=xCorr;
			this.spectrumIndex=spectrumIndex;
			this.numberOfPeaks=numberOfPeaks;
			this.sumCorrelation=sumCorrelation;
			this.foundLocalizingIons=foundLocalizingIons;
		}
	}

	/**
	 * extracts chromatograms for all peaks at once
	 * @author searleb
	 *
	 */
	private class ChromatogramMap {
		float[] retentionTimes;
		TDoubleObjectHashMap<float[]> chromatogramsByTargetMass;
		public ChromatogramMap(double[] targetMasses, int peakIndex) {
			List<Stripe> localStripes=stripes.subList(Math.max(0, peakIndex-movingAverageLength), Math.min(stripes.size(), peakIndex+movingAverageLength));
			TFloatArrayList[] chromatograms=new TFloatArrayList[targetMasses.length];
			for (int i = 0; i < chromatograms.length; i++) {
				chromatograms[i]=new TFloatArrayList();
			}
			
			TFloatArrayList rtList=new TFloatArrayList();
			for (Stripe spectrum : localStripes) {
				rtList.add(spectrum.getScanStartTime());
				float[] integratedIntensities=parameters.getFragmentTolerance().getIntegratedIntensities(spectrum.getMassArray(), spectrum.getIntensityArray(), targetMasses);
				for (int i = 0; i < chromatograms.length; i++) {
					chromatograms[i].add(integratedIntensities[i]);
				}
			}
			retentionTimes=rtList.toArray();

			chromatogramsByTargetMass=new TDoubleObjectHashMap<>();
			for (int i = 0; i < chromatograms.length; i++) {
				chromatogramsByTargetMass.put(targetMasses[i], chromatograms[i].toArray());
			}
		}
		
		/**
		 * assumes that target masses are bitwise identical doubles!
		 * @param target
		 * @return
		 */
		public float[] getChromatogram(double target) {
			return chromatogramsByTargetMass.get(target);
		}
		
		public ArrayList<float[]> getChromatograms(double[] targets) {
			ArrayList<float[]> chromatograms=new ArrayList<>();
			for (int i = 0; i < targets.length; i++) {
				chromatograms.add(getChromatogram(targets[i]));
			}
			return chromatograms;
		}
		
		public float[] getRetentionTimes() {
			return retentionTimes;
		}
	}

	private static XCorrLibraryEntry getXCorrEntry(LibraryEntry entry, SearchParameters parameters) {
		XCorrLibraryEntry xcordiaEntry;
		if (entry instanceof XCorrLibraryEntry) {
			xcordiaEntry=(XCorrLibraryEntry)entry;
		} else {
			FastaPeptideEntry peptide=new FastaPeptideEntry(entry.getSource(), entry.getAccessions(), entry.getPeptideModSeq());
			xcordiaEntry=XCorrLibraryEntry.generateEntry(false, peptide, entry.getPrecursorCharge(), parameters);
		}
		xcordiaEntry.init();
		return xcordiaEntry;
	}
}