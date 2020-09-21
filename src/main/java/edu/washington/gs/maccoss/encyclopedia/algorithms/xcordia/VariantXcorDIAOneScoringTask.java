package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
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
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.FragmentIonBlacklist;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoPermuter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusOneScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.ThesaurusOneScoringTask.LocalizableForm;
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
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
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
	private final PeptideModification localizingModification;
	
	public VariantXcorDIAOneScoringTask(PSMScorer scorer, BackgroundFrequencyInterface background, ArrayList<LibraryEntry> entries, 
			ArrayList<FragmentScan> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<PeptideScoringResult> resultsQueue,
			BlockingQueue<ModificationLocalizationData> localizationQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.background=background;
		this.dutyCycle=dutyCycle;
		this.localizationQueue=localizationQueue;
		this.minimumScore=-Log.log10(parameters.getPercolatorThreshold());
		this.movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		if (parameters.getLocalizingModification().isPresent()) {
			this.localizingModification=parameters.getLocalizingModification().get();
		} else {
			this.localizingModification=PeptideModification.polymorphism; 
		}
		
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
	
	private void processPeptide(ArrayList<LibraryEntry> seedEntries) {
		LibraryEntry firstEntry=seedEntries.get(0);
		ArrayList<String> peptideModSeqs=PhosphoPermuter.getPermutations(firstEntry.getPeptideModSeq(), localizingModification, parameters.getAAConstants());

		byte precursorCharge=firstEntry.getPrecursorCharge();

		HashMap<String, double[]> ionsByPeptide=new HashMap<>();
		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
		for (String peptideModSeq : peptideModSeqs) {
			FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants());
			entryMap.put(peptideModSeq, model);
			ionsByPeptide.put(peptideModSeq, model.getPrimaryIons(parameters.getFragType(), firstEntry.getPrecursorCharge(), true));
		}

		HashMap<String, Range> localizedIsoforms=new HashMap<>();
		HashSet<String> unlocalizedIsoforms=new HashSet<>();
		HashMap<String, LocalizableForm> allIsoforms=new HashMap<>();
		HashMap<String, float[]> predictedIsotopeDistributions=new HashMap<>();
		for (String peptideModSeq : peptideModSeqs) {
			LocalizableForm form=ThesaurusOneScoringTask.getLocalizedForm(peptideModSeq, firstEntry.getPrecursorCharge(), entryMap, ionsByPeptide, seedEntries, parameters);
			if (form!=null) {
				allIsoforms.put(peptideModSeq, form);
				unlocalizedIsoforms.add(peptideModSeq);
			}

			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(peptideModSeq, parameters.getAAConstants());
			predictedIsotopeDistributions.put(peptideModSeq, predictedIsotopeDistribution);
		}

		@SuppressWarnings("unchecked")
		ArrayList<Spectrum>[] globalScans=new ArrayList[] {PhosphoLocalizer.getScanSubsetFromStripes(-Float.MAX_VALUE, Float.MAX_VALUE, super.stripes)};
		@SuppressWarnings("unchecked")
		HashMap<String, Float[]>[] globalPrimaryScores=new HashMap[globalScans.length]; // scores are re-used unless they fall in the RT range of a previously localized form
		for (int i=0; i<globalPrimaryScores.length; i++) {
			globalPrimaryScores[i]=new HashMap<>();
		}
		@SuppressWarnings("unchecked")
		HashMap<String, boolean[]>[] globalAlreadyConsideredScores=new HashMap[globalScans.length];
		for (int i=0; i<globalAlreadyConsideredScores.length; i++) {
			globalAlreadyConsideredScores[i]=new HashMap<>();
		}
		
		HashMap<String, Range> blacklistedScanRanges=new HashMap<>();
		
		PeptideScoringResult bestNonlocalizedResult=null;
		ModificationLocalizationData bestNonlocalizedData=null;
		boolean anyLocalized=false;
		
		HashMap<String, ScoredIndex> bestIndicies=new HashMap<>();
		
		while (unlocalizedIsoforms.size()>0) {
			Collection<Range> blacklistedScans=blacklistedScanRanges.values();
			String bestPeptideModSeq=null;
			ScoredIndex bestIndex=null;
			LocalizableForm bestForm=null;
			ArrayList<Spectrum> scans=null;
			HashMap<String, Float[]> primaryScores=null; // scores are re-used unless they fall in the RT range of a previously localized form
			HashMap<String, boolean[]> alreadyConsideredScores=null;
			
			// get the highest scoring RT point for all peptide sequences
			for (Entry<String, LocalizableForm> entry : allIsoforms.entrySet()) {
				String peptideModSeq=entry.getKey();
				LocalizableForm localizedForm=entry.getValue();
				
				XCorrLibraryEntry localizedEntry=(XCorrLibraryEntry)localizedForm.localizedEntry;

				for (int peak=0; peak<globalScans.length; peak++) {
					Float[] primary=globalPrimaryScores[peak].get(peptideModSeq);
					boolean[] alreadyConsidered=globalAlreadyConsideredScores[peak].get(peptideModSeq);

					if (primary==null) {
						primary=new Float[globalScans[peak].size()];
						alreadyConsidered=new boolean[primary.length];
						globalPrimaryScores[peak].put(peptideModSeq, primary);
						globalAlreadyConsideredScores[peak].put(peptideModSeq, alreadyConsidered);
					}
					ScoredIndex score=updateScores(globalScans[peak], localizedEntry, predictedIsotopeDistributions.get(peptideModSeq), primary, alreadyConsidered, blacklistedScans, parameters);
					
					ScoredIndex previousBest=bestIndicies.get(peptideModSeq);
					if (previousBest==null) {
						previousBest=score;
						bestIndicies.put(peptideModSeq, score);
					}

					if (score!=null&&unlocalizedIsoforms.contains(peptideModSeq)&&(bestIndex==null||score.x>bestIndex.x)&&(previousBest.x*ThesaurusOneScoringTask.MAXIMUM_DELTA_FROM_BEST_SCORE<score.x)) {
						bestIndex=score;
						bestPeptideModSeq=peptideModSeq;
						bestForm=localizedForm;
						
						scans=globalScans[peak];
						primaryScores=globalPrimaryScores[peak];
						alreadyConsideredScores=globalAlreadyConsideredScores[peak];
					}
				}
			}
			//System.out.println("CHECK: "+bestPeptideModSeq+"\t"+bestIndex.x+"\t"+(scans.get(bestIndex.y).getScanStartTime()/60f)+" (total scans: "+scans.size()+"), attempts:"+numberOfAttempts.get(bestPeptideModSeq)); //FIXME
			if (bestIndex==null) {
				// no more good matches
				break;
			}
			ScoredIndex previousScoreIndex=bestIndicies.get(bestPeptideModSeq);
			if (previousScoreIndex!=null&&previousScoreIndex.x*ThesaurusOneScoringTask.MAXIMUM_DELTA_FROM_BEST_SCORE>=bestIndex.x) {
				unlocalizedIsoforms.remove(bestPeptideModSeq); // bad score, not worth considering
				continue;
			}
			
//			if (first) { //FIXME
//				ArrayList<XYTraceInterface> traces=new ArrayList<>();
//				for (int peak=0; peak<globalScans.length; peak++) {
//					float[] rts=new float[scans.size()];
//					for (int i = 0; i < rts.length; i++) {
//						rts[i]=globalScans[peak].get(i).getScanStartTime()/60f;
//					}
//					for (Entry<String, Float[]> scores : globalPrimaryScores[peak].entrySet()) {
//						traces.add(new XYTrace(rts, General.toFloatArray(scores.getValue()), GraphType.boldline, seedEntries.get(peak).getPeptideModSeq()+"/"+scores.getKey()));
//					}
//				}
//				Charter.launchChart("RT", "Score", true, traces.toArray(new XYTraceInterface[traces.size()]));
//				first=false;
//			}
			
			// check localization ions versus that sequence
			float apexRT=scans.get(bestIndex.y).getScanStartTime();
			// use stripes here in case we're on the border
			Range generalPeakArea=new Range(apexRT-parameters.getExpectedPeakWidth(), apexRT+parameters.getExpectedPeakWidth());
			
			ArrayList<Spectrum> stripeSubset=PhosphoLocalizer.getScanSubsetFromStripes(generalPeakArea.getStart(), generalPeakArea.getStop(), stripes); 

			// get next best match at that RT
			FragmentIon[] bestLocalizingIons=null;
			Pair<FragmentScan, Float> bestLocalizedStripe=null;
			
			for (Entry<String, Float[]> entry : primaryScores.entrySet()) {
				String peptideModSeq=entry.getKey();
				if (peptideModSeq!=bestPeptideModSeq) {
					float score=entry.getValue()[bestIndex.y];
					Range previousLocalization=localizedIsoforms.get(peptideModSeq);
					// if we've previously localized this peptide 
					// AND it wasn't at this location 
					// AND it scores WORSE than the current peptide (i.e. it's not a peptide that's eluting twice)
					// THEN (and only then) we can confidently say it isn't here, so we can ignore it
					if (previousLocalization!=null&&!previousLocalization.contains(apexRT)&&score<bestIndex.x) {
						continue;
					}
					
					FragmentIon[] localizingIons=ThesaurusOneScoringTask.getUniqueFragmentIons(entryMap.get(bestPeptideModSeq), entryMap.get(peptideModSeq), precursorCharge, parameters);
					
					Pair<FragmentScan, Float> localizedStripe = ThesaurusOneScoringTask.getBestLocalizationStripe(parameters, dutyCycle, background, bestForm.localizedEntry, localizingIons, stripeSubset);
					//System.out.println("Testing "+bestPeptideModSeq+" ("+bestIndex.x+") vs "+peptideModSeq+" ("+score+"): localization: "+localizedStripe.y);
					if (bestLocalizedStripe==null||bestLocalizedStripe.y>localizedStripe.y) {
						// keep the lowest localization scoring form! (the closest to the form we're considering)
						bestLocalizedStripe=localizedStripe;
						bestLocalizingIons=localizingIons;
					}
				}
			}	
			
			if (bestLocalizedStripe==null) {
				bestLocalizingIons=bestForm.allIons;
				bestLocalizedStripe = ThesaurusOneScoringTask.getBestLocalizationStripe(parameters, dutyCycle, background, bestForm.localizedEntry, bestLocalizingIons, stripeSubset);
			}
			
			AmbiguousPeptideModSeq ambiPeptideModSeq=AmbiguousPeptideModSeq.getUnambigous(bestPeptideModSeq, localizingModification, parameters.getAAConstants(), "");
			Triplet<ModificationLocalizationData, FragmentScan, Range> locData=ThesaurusOneScoringTask.generateLocalizationData(false, minimumScore, parameters, localizingModification, bestForm.localizedEntry,
					ambiPeptideModSeq, bestLocalizingIons, bestForm.allIons, takenIdentifiedIons, stripeSubset, bestLocalizedStripe);
			if (locData==null) continue;
			
			// if localized, then keep and remove from localizedForms
			ModificationLocalizationData data=locData.x;
			FragmentScan apex=locData.y;
			Range peakRange=locData.z.addBuffer(dutyCycle);

			XCorrStripe xcordiaStripe;
			if (apex instanceof XCorrStripe) {
				xcordiaStripe=(XCorrStripe)apex;
			} else {
				xcordiaStripe=new XCorrStripe(apex, parameters);
			}
			float[] predictedIsotopeDistribution=predictedIsotopeDistributions.get(bestForm.localizedEntry.getPeptideModSeq());
			float score=scorer.score(bestForm.localizedEntry, xcordiaStripe, predictedIsotopeDistribution, precursors);
			boolean replaceBestNonLocalizedResult = bestNonlocalizedResult==null||(score>=bestNonlocalizedResult.getBestScore()&&bestNonlocalizedData.getLocalizationScore()<data.getLocalizationScore());
			//System.out.println(bestPeptideModSeq+" --> "+apex.getScanStartTime()/60f+" min, loc:"+locData.x.getLocalizationScore()+" vs "+bestLocalizedStripe.y+", score:"+score+" vs "+(bestNonlocalizedResult==null?"null":Float.toString(bestNonlocalizedResult.getBestScore()))+" ("+replaceBestNonLocalizedResult+")");
			if (replaceBestNonLocalizedResult||data.isLocalized()) {
				TFloatFloatHashMap scoreByRTMap=new TFloatFloatHashMap();
				Float[] primaryScoreArray=primaryScores.get(bestPeptideModSeq);
				for (int i=0; i<scans.size(); i++) {
					scoreByRTMap.put(scans.get(i).getScanStartTime(), primaryScoreArray[i]);
				}
				EValueCalculator calculator=new EValueCalculator(scoreByRTMap, 0.1f, 0.1f);
				float[] auxScoreArray=scorer.auxScore(bestForm.localizedEntry, apex, predictedIsotopeDistribution, precursors);
	
				float evalue=calculator.getNegLnEValue(score);
				if (Float.isNaN(evalue)) {
					evalue=-1.0f;
				}

				PeptideScoringResult result=new PeptideScoringResult(bestForm.localizedEntry);
				result.addStripe(score, General.concatenate(auxScoreArray, evalue, data.getLocalizationScore(), data.getNumIdentificationPeaks()), apex);

				if (replaceBestNonLocalizedResult) {
					bestNonlocalizedResult=result;
					bestNonlocalizedData=data;
				}
				
				if (data.isLocalized()) {
					anyLocalized=true;
					resultsQueue.add(result);

					if (!bestForm.localizedEntry.isDecoy()) {
						// don't bother logging decoys
						localizationQueue.add(data);
					}
				}
			}

			if (data.isLocalized()) {
				// we've found this peptide so we don't need to keep this section blacklisted
				blacklistedScanRanges.remove(bestPeptideModSeq);
			} else {
				if (!blacklistedScanRanges.containsKey(bestPeptideModSeq)) {
					// if we've already got a better RT peak blacklisted then don't bother with this one
					blacklistedScanRanges.put(bestPeptideModSeq, generalPeakArea);
				}				
			}

			//System.out.println("Blocking off "+(peakRange.getStart()/60f)+" to "+(peakRange.getStop()/60f)+" for "+bestPeptideModSeq+" --> "+data.isLocalized()+", "+data.getLocalizationScore());
			for (Ion target : bestLocalizingIons) {
				takenIdentifiedIons.addIonToBlacklist(target.getMass(), peakRange);
			}
			// null out scores from taken ions
			for (int i=0; i<scans.size(); i++) {
				if (peakRange.contains(scans.get(i).getScanStartTime())) {
					for (Entry<String, Float[]> entry : primaryScores.entrySet()) {
						String peptideModSeq = entry.getKey();
						if (peptideModSeq!=bestPeptideModSeq) {
							if (unlocalizedIsoforms.contains(peptideModSeq)) {// &&entry.getValue()[i]!=THIS_PEPTIDE_IS_NOT_HERE) {
								entry.getValue()[i]=null;
							}
						} else if (!data.isLocalized()) {
							// if not localized and around the peak shape area
							alreadyConsideredScores.get(peptideModSeq)[i]=true;
						}
					}
				} else {
					if (!data.isLocalized()&&generalPeakArea.contains(scans.get(i).getScanStartTime())) {
						// if not localized and anywhere near the apex
						alreadyConsideredScores.get(bestPeptideModSeq)[i]=true;
					}
				}
			}
			
			if (data.isLocalized()) {
				localizedIsoforms.put(bestPeptideModSeq, peakRange);
			}

			// try only one time!
			unlocalizedIsoforms.remove(bestPeptideModSeq); // should we only do this if we can actually localize the peak?
		}
		
		if (!anyLocalized&&bestNonlocalizedResult!=null) {
			resultsQueue.add(bestNonlocalizedResult);

			if (!bestNonlocalizedResult.getEntry().isDecoy()) {
				// don't bother logging decoys
				localizationQueue.add(bestNonlocalizedData);
			}
		}
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
		
		PeptideScoringResult result=new PeptideScoringResult(timepoint.xcordiaEntry);
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
	


	private ScoredIndex updateScores(ArrayList<Spectrum> scans, XCorrLibraryEntry localizedEntry, float[] predictedIsotopeDistribution, Float[] primary, boolean[] alreadyConsidered, Collection<Range> blacklistedRanges, SearchParameters parameters) {
		float bestScore=-Float.MAX_VALUE;
		int bestIndex=0;
		for (int i=0; i<scans.size(); i++) {
			Spectrum stripe=scans.get(i);
			XCorrStripe xcordiaStripe;
			if (stripe instanceof XCorrStripe) {
				xcordiaStripe=(XCorrStripe)stripe;
			} else {
				xcordiaStripe=new XCorrStripe((FragmentScan)stripe, parameters);
			}
			
			if (primary[i]==null) { 
				primary[i]=scorer.score(localizedEntry, xcordiaStripe, predictedIsotopeDistribution, precursors);
			}
			if ((!alreadyConsidered[i])&&(bestScore<primary[i])) {
				boolean alreadyTaken=false;
				for (Range range : blacklistedRanges) {
					if (range.contains(stripe.getScanStartTime())) {
						alreadyTaken=true;
					}
				}
				if (!alreadyTaken) {
					bestIndex=i;
					bestScore=primary[i];
				}
			}
		}
		
		return new ScoredIndex(bestScore, bestIndex);
	}
}
