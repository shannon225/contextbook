package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.FragmentIonBlacklist;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific.VariantFastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class VariantXcorDIAOneScoringTask extends AbstractLibraryScoringTask {
	private final float dutyCycle;
	private final BackgroundFrequencyInterface background;
	private final BlockingQueue<ModificationLocalizationData> localizationQueue;
	private final float minimumScore;
	private final int movingAverageLength;
	
	public VariantXcorDIAOneScoringTask(PSMScorer scorer, BackgroundFrequencyInterface background, ArrayList<LibraryEntry> entries, 
			ArrayList<FragmentScan> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue,
			BlockingQueue<ModificationLocalizationData> localizationQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.background=background;
		this.dutyCycle=dutyCycle;
		this.localizationQueue=localizationQueue;
		this.minimumScore=-Log.log10(parameters.getPercolatorThreshold());
		this.movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		
		assert(scorer instanceof XCorDIAOneScorer);
	}

	@Override
	protected Nothing process() {
		// separate targets from decoys and process in batches
		ArrayList<LibraryEntry> targetBatch=new ArrayList<>();
		ArrayList<LibraryEntry> decoyBatch=new ArrayList<>();
		for (LibraryEntry entry : super.entries) {
			XCorrLibraryEntry xcordiaEntry = getXCorrEntry(entry);
			if (entry.isDecoy()) {
				decoyBatch.add(xcordiaEntry);
			} else {
				targetBatch.add(xcordiaEntry);
			}
		}
		processBatch(targetBatch);
		processBatch(decoyBatch);

		return Nothing.NOTHING;
	}

	// NOTE: assumes all peptides in entries are related and that decoys are run separately than targets
	private void processBatch(ArrayList<LibraryEntry> entryBatch) {
		HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
		for (LibraryEntry entry : entryBatch) {
			FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
			modelBatch.put(entry.getPeptideModSeq(), model);
		}
		
		ScoredTimepoint bestTimepoint=null;
		boolean reportedScores=false;
		
		for (LibraryEntry entry : entryBatch) {
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(entry.getPeptideModSeq(), parameters.getAAConstants());

			// determine index of best starting point 
			XCorrLibraryEntry xcordiaEntry = getXCorrEntry(entry);
			float[] primary = scoreEntryAcrossTime(xcordiaEntry, predictedIsotopeDistribution);
			
			TFloatFloatHashMap map=new TFloatFloatHashMap();
			for (int i=0; i<primary.length; i++) {
				map.put(i, primary[i]); // rt=index
			}
			float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);
			
			EValueCalculator calculator=new EValueCalculator(map, 0.1f, 0.1f);
			int index=Math.round(calculator.getMaxRT()); // rt=index
			
			IntRange indexRange = getPeakRange(index);
			
			// get localizing ions			
			FragmentIon[] targets=PhosphoLocalizer.getUniqueFragmentIons(entry.getPeptideModSeq(), entry.getPrecursorCharge(), modelBatch, parameters);
			double[] ions=FragmentIon.getMasses(targets);
			float[] frequencies=background.getFrequencies(ions, entry.getPrecursorMZ(), parameters.getFragmentTolerance());
			
			// run localization scoring on target ions
			int bestLocalizationIndex=0;
			float maxLocalizationScore=-Float.MAX_VALUE;
			float maxLocalizedPrimary=-Float.MAX_VALUE;
			float[] negLogProbsSiteSpecific=new float[indexRange.getLength()];
			for (int i = indexRange.getStart(); i <= indexRange.getStop(); i++) {
				float rawScore = PhosphoLocalizer.score(parameters, ions, targets, frequencies, stripes.get(i), true);
				negLogProbsSiteSpecific[i-indexRange.getStart()]=rawScore;
				if (rawScore>maxLocalizationScore||(rawScore==maxLocalizationScore&&averagePrimary[i]>maxLocalizedPrimary)) {
					maxLocalizationScore=rawScore;
					maxLocalizedPrimary=averagePrimary[i];
					bestLocalizationIndex=i;
				}
			}
			negLogProbsSiteSpecific=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, Math.round(parameters.getExpectedPeakWidth()/(dutyCycle)));//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);
			
			// calculate final scoring on best localization index 
			ScoredIndex scanIndex=new ScoredIndex(primary[bestLocalizationIndex], bestLocalizationIndex);
			FragmentIon[] allIons=modelBatch.get(entry.getPeptideModSeq()).getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true);

			ArrayList<Spectrum> stripeSubset=PhosphoLocalizer.getScanSubsetFromStripes(stripes.get(indexRange.getStart()).getScanStartTime(), stripes.get(indexRange.getStop()).getScanStartTime()+parameters.getExpectedPeakWidth(), stripes);
			Pair<ModificationLocalizationData, Integer> pair=getLocalizationData(xcordiaEntry.getPeptide(), stripes.get(bestLocalizationIndex), entry.getPeptideModSeq(), entry.getPrecursorCharge(), 
					minimumScore, maxLocalizationScore, targets, allIons, stripeSubset);
			
			localizationQueue.add(pair.x);
			ScoredTimepoint timepoint = new ScoredTimepoint(xcordiaEntry, predictedIsotopeDistribution, calculator, scanIndex, entryBatch.size()>1, pair.y);
			
			if (bestTimepoint==null||timepoint.scoredIndex.x>bestTimepoint.scoredIndex.x) {
				bestTimepoint=timepoint;
			}
			if (pair.x.isLocalized()) {
				// report every localized timepoint
				finalScoreTimepoint(timepoint);
				reportedScores=true;
			}
		}
		
		// always report something
		if (!reportedScores&&bestTimepoint!=null) {
			finalScoreTimepoint(bestTimepoint);
		}
	}

	private IntRange getPeakRange(int index) {
		// find region around best scoring scan
		int startIndex=index;
		int stopIndex=index;
		for (int i = index-1; i>=0; i--) {
			if (stripes.get(index).getScanStartTime()-parameters.getExpectedPeakWidth()<stripes.get(i).getScanStartTime()) {
				startIndex=i;
			} else {
				break;
			}
		}
		for (int i = index+1; i < stripes.size(); i++) {
			if (stripes.get(index).getScanStartTime()+parameters.getExpectedPeakWidth()>stripes.get(i).getScanStartTime()) {
				stopIndex=i;
			} else {
				break;
			}
		}
		IntRange indexRange=new IntRange(startIndex, stopIndex);
		return indexRange;
	}
	
	FragmentIonBlacklist takenIdentifiedIons=new FragmentIonBlacklist(parameters.getFragmentTolerance()); // not necessary
	private Pair<ModificationLocalizationData, Integer> getLocalizationData(FastaPeptideEntry peptide, FragmentScan apex, String peptideModSeq, byte precursorCharge, float minimumScore, float bestLocalizationScore, FragmentIon[] targetIons, FragmentIon[] allIons, ArrayList<Spectrum> stripeSubset) {
		int targetNumFragments=Math.max(parameters.getMinNumOfQuantitativePeaks(), 3);
		
		//Range peakRange=new Range(apex.getScanStartTime(), apex.getScanStartTime());
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" --> ("+targetIons.length+") "+bestLocalizationScore+" localization score at "+(apex.getScanStartTime()/60f)+" min"); //FIXME 

		boolean isLocalized=false;
		boolean isSiteSpecific=false;
		boolean isCompletelyAmbiguous=false;
		float localizationIntensity=0.0f;
		float totalIntensity=0.0f;
		int numIdentificationPeaks=0;
		int numberOfMods=0;
		float apexRT=apex.getScanStartTime();
		ArrayList<FragmentIon> wellShapedIons=new ArrayList<FragmentIon>();
		
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score)\tNOT GOOD ENOUGH");
		if (bestLocalizationScore>=minimumScore) {
			// generate quant data from localizing ions only
			TransitionRefinementData quantData=PhosphoLocalizer.quantifyPeptide(parameters, peptideModSeq, precursorCharge, targetIons, apex.getScanStartTime(),
					stripeSubset, takenIdentifiedIons, Optional.ofNullable((float[]) null));
			apexRT=quantData.getApexRT();
			//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score)\ta)"+(quantData!=null));
			if (quantData!=null) {
				//peakRange=quantData.getRange();
				float[] intensities=quantData.getIntegrationArray();
				float[] correlations=quantData.getCorrelationArray();
				FragmentIon[] consideredIons=(FragmentIon[])quantData.getFragmentMassArray();
				for (int i=0; i<consideredIons.length; i++) {
					wellShapedIons.add(consideredIons[i]);
					localizationIntensity+=intensities[i];
				}

				// calculate quant data for all ions
				float[] medianChromatogram=quantData.getMedianChromatogram();
				TransitionRefinementData allQuantData=PhosphoLocalizer.quantifyPeptide(parameters, peptideModSeq, precursorCharge, allIons, apex.getScanStartTime(),
						stripeSubset, takenIdentifiedIons, Optional.of(medianChromatogram));
				//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score)\tb)"+(allQuantData!=null));
				if (allQuantData!=null) {
					//peakRange=allQuantData.getRange();
					intensities=allQuantData.getIntegrationArray();
					correlations=allQuantData.getCorrelationArray();
					for (int i=0; i<correlations.length; i++) {
						if (correlations[i]>=TransitionRefiner.identificationCorrelationThreshold) {
							numIdentificationPeaks++;
							totalIntensity+=intensities[i];
						}
					}
					isCompletelyAmbiguous=false;
					isLocalized=wellShapedIons.size()>0&&numIdentificationPeaks>=targetNumFragments&&!isCompletelyAmbiguous;
					isSiteSpecific=isLocalized;
					
					//System.out.println("\tLocalized "+isSiteSpecific+" for "+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score, "+numIdentificationPeaks+"/"+correlations.length+" peaks)"); // FIXME
				}
			}
		}
		String annotation;
		if (peptide instanceof VariantFastaPeptideEntry) {
			annotation=((VariantFastaPeptideEntry)peptide).getVariant().toString();
		} else {
			annotation="canonical";
		}
		AmbiguousPeptideModSeq ambiguousPeptideModSeq=AmbiguousPeptideModSeq.getUnambigous(peptideModSeq, PeptideModification.polymorphism, parameters.getAAConstants(), annotation);
		
		ModificationLocalizationData modData=new ModificationLocalizationData(ambiguousPeptideModSeq, apexRT, bestLocalizationScore, numIdentificationPeaks, numberOfMods, isSiteSpecific, isLocalized, isCompletelyAmbiguous, wellShapedIons.toArray(new FragmentIon[wellShapedIons.size()]), localizationIntensity, totalIntensity);
		return new Pair<ModificationLocalizationData, Integer>(modData, numIdentificationPeaks);
	}

	private XCorrLibraryEntry getXCorrEntry(LibraryEntry entry) {
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
	
	private class ScoredTimepoint {
		XCorrLibraryEntry xcordiaEntry;
		float[] predictedIsotopeDistribution;
		EValueCalculator calculator;
		ScoredIndex scoredIndex;
		boolean neededToLocalize;
		int numberOfWellShapedIons;
		public ScoredTimepoint(XCorrLibraryEntry xcordiaEntry, float[] predictedIsotopeDistribution,
				EValueCalculator calculator, ScoredIndex scoredIndex, boolean neededToLocalize,
				int numberOfWellShapedIons) {
			this.xcordiaEntry = xcordiaEntry;
			this.predictedIsotopeDistribution = predictedIsotopeDistribution;
			this.calculator = calculator;
			this.scoredIndex = scoredIndex;
			this.neededToLocalize = neededToLocalize;
			this.numberOfWellShapedIons = numberOfWellShapedIons;
		}
	}

	private void finalScoreTimepoint(ScoredTimepoint timepoint) {
		FragmentScan stripe=super.stripes.get(timepoint.scoredIndex.y);
		float[] auxScoreArray=scorer.auxScore(timepoint.xcordiaEntry, stripe, timepoint.predictedIsotopeDistribution, precursors);
		float evalue=timepoint.calculator.getNegLnEValue(timepoint.scoredIndex.x);
		if (Float.isNaN(evalue)) {
			evalue=-1.0f;
		}
		
		AbstractScoringResult result=new PeptideScoringResult(timepoint.xcordiaEntry);
		result.addStripe(timepoint.scoredIndex.x, General.concatenate(auxScoreArray, timepoint.scoredIndex.x, evalue, timepoint.neededToLocalize?1:0, timepoint.numberOfWellShapedIons), stripe);
		resultsQueue.add(result);
	}

	private float[] scoreEntryAcrossTime(XCorrLibraryEntry xcordiaEntry, float[] predictedIsotopeDistribution) {
		float[] primary=new float[super.stripes.size()];
		float[] rts=new float[super.stripes.size()];
		for (int i=0; i<super.stripes.size(); i++) {
			FragmentScan stripe=super.stripes.get(i);
			XCorrStripe xcordiaStripe;
			if (stripe instanceof XCorrStripe) {
				xcordiaStripe=(XCorrStripe)stripe;
			} else {
				xcordiaStripe=new XCorrStripe(stripe, parameters);
			}
			primary[i]=scorer.score(xcordiaEntry, xcordiaStripe, predictedIsotopeDistribution, precursors);
			
			rts[i]=xcordiaStripe.getScanStartTime();
		}
		return primary;
	}
}
