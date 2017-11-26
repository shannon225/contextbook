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
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.BackgroundFrequencyInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.FragmentIonBlacklist;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
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
			ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue,
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
			if (entry.isDecoy()) {
				decoyBatch.add(entry);
			} else {
				targetBatch.add(entry);
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
			FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
			modelBatch.put(entry.getPeptideModSeq(), model);
		}
		
		for (LibraryEntry entry : entryBatch) {
			XCorrLibraryEntry xcordiaEntry = getXCorrEntry(entry);
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(xcordiaEntry.getPeptideModSeq(), parameters.getAAConstants());

			// determine index of best starting point 
			float[] primary = scoreEntryAcrossTime(xcordiaEntry, predictedIsotopeDistribution);
			TFloatFloatHashMap map=new TFloatFloatHashMap();
			for (int i=0; i<primary.length; i++) {
				map.put(i, primary[i]); // rt=index
			}
			float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);
			
			EValueCalculator calculator=new EValueCalculator(map);
			int index=Math.round(calculator.getMaxRT()); // rt=index

			// no need to localize
			if (entryBatch.size()==1) {
				ScoredIndex scanIndex=new ScoredIndex(primary[index], index);
				finalScoreTimepoint(xcordiaEntry, predictedIsotopeDistribution, calculator, scanIndex);
				return;
			}
			IntRange indexRange = getPeakRange(index);
			
			// get localizing ions			
			FragmentIon[] targets=PhosphoLocalizer.getUniqueFragmentIons(entry.getPeptideModSeq(), entry.getPrecursorCharge(), modelBatch, parameters);
			double[] ions=FragmentIon.getMasses(targets);
			float[] frequencies=background.getFrequencies(ions, entry.getPrecursorMZ(), parameters.getFragmentTolerance());
			
			// run localization scoring on target ions
			int bestIndex=0;
			float maxRawScore=-Float.MAX_VALUE;
			float maxRawPrimary=-Float.MAX_VALUE;
			float[] negLogProbsSiteSpecific=new float[indexRange.getLength()];
			for (int i = indexRange.getStart(); i <= indexRange.getStop(); i++) {
				float rawScore = PhosphoLocalizer.score(parameters, ions, targets, frequencies, stripes.get(i), true);
				negLogProbsSiteSpecific[i-indexRange.getStart()]=rawScore;
				if (rawScore>maxRawScore||(rawScore==maxRawScore&&averagePrimary[i]>maxRawPrimary)) {
					maxRawScore=rawScore;
					maxRawPrimary=averagePrimary[i];
					bestIndex=i;
				}
			}
			negLogProbsSiteSpecific=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, Math.round(parameters.getExpectedPeakWidth()/(dutyCycle)));//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);
			
			// calculate final scoring on best localization index 
			ScoredIndex scanIndex=new ScoredIndex(primary[bestIndex], bestIndex);
			finalScoreTimepoint(xcordiaEntry, predictedIsotopeDistribution, calculator, scanIndex);
			if (!xcordiaEntry.isDecoy()) {
				// don't bother logging decoys
				FragmentIon[] allIons=modelBatch.get(xcordiaEntry.getPeptideModSeq()).getPrimaryIonObjects(parameters.getFragType(), xcordiaEntry.getPrecursorCharge());

				ArrayList<Spectrum> stripeSubset=PhosphoLocalizer.getScanSubsetFromStripes(stripes.get(indexRange.getStart()).getScanStartTime(), stripes.get(indexRange.getStop()).getScanStartTime()+parameters.getExpectedPeakWidth(), stripes);
				ModificationLocalizationData data=getLocalizationData(stripes.get(bestIndex), xcordiaEntry.getPeptideModSeq(), xcordiaEntry.getPrecursorCharge(), 
						minimumScore, maxRawScore, targets, allIons, stripeSubset);
				
				localizationQueue.add(data);
			}
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
	private ModificationLocalizationData getLocalizationData(Stripe apex, String peptideModSeq, byte precursorCharge, float minimumScore, float bestLocalizationScore, FragmentIon[] targetIons, FragmentIon[] allIons, ArrayList<Spectrum> stripeSubset) {
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
				FragmentIon[] consideredIons=quantData.getFragmentMassArray();
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
		AmbiguousPeptideModSeq ambiguousPeptideModSeq=AmbiguousPeptideModSeq.getUnambigous(peptideModSeq, PeptideModification.polymorphism, parameters.getAAConstants());
		
		ModificationLocalizationData modData=new ModificationLocalizationData(ambiguousPeptideModSeq, apexRT, bestLocalizationScore, numberOfMods, isSiteSpecific, isLocalized, isCompletelyAmbiguous, wellShapedIons.toArray(new FragmentIon[wellShapedIons.size()]), localizationIntensity, totalIntensity);
		return modData;
	}

	private XCorrLibraryEntry getXCorrEntry(LibraryEntry entry) {
		XCorrLibraryEntry xcordiaEntry;
		if (entry instanceof XCorrLibraryEntry) {
			xcordiaEntry=(XCorrLibraryEntry)entry;
		} else {
			xcordiaEntry=XCorrLibraryEntry.generateEntry(false, entry.getSource(), entry.getAccessions(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), parameters);
		}
		xcordiaEntry.init();
		return xcordiaEntry;
	}

	private void finalScoreTimepoint(XCorrLibraryEntry xcordiaEntry, float[] predictedIsotopeDistribution,
			EValueCalculator calculator, ScoredIndex scoredIndex) {
		Stripe stripe=super.stripes.get(scoredIndex.y);
		float[] auxScoreArray=scorer.auxScore(xcordiaEntry, stripe, predictedIsotopeDistribution, precursors);
		float evalue=calculator.getNegLog10EValue(scoredIndex.x);
		if (Float.isNaN(evalue)) {
			evalue=-1.0f;
		}
		
		PeptideScoringResult result=new PeptideScoringResult(xcordiaEntry);
		result.addStripe(scoredIndex.x, General.concatenate(auxScoreArray, scoredIndex.x, evalue), stripe);
		resultsQueue.add(result);
	}

	private float[] scoreEntryAcrossTime(XCorrLibraryEntry xcordiaEntry, float[] predictedIsotopeDistribution) {
		float[] primary=new float[super.stripes.size()];
		float[] rts=new float[super.stripes.size()];
		for (int i=0; i<super.stripes.size(); i++) {
			Stripe stripe=super.stripes.get(i);
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
