package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AuxillaryPSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class CASiLOneScoringTask extends AbstractLibraryScoringTask {
	private static final int peaksKept=-1; // only keep the first peak
	
	private final PhosphoLocalizer localizer;
	private final float dutyCycle;
	private final ScoringBreadthType breadth;
	private final PeptideModification localizingModification;
	private final BlockingQueue<ModificationLocalizationData> localizationQueue;
	private final float minimumScore;
	
	public CASiLOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, 
			PhosphoLocalizer localizer, BlockingQueue<PeptideScoringResult> resultsQueue, BlockingQueue<ModificationLocalizationData> localizationQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.dutyCycle=dutyCycle;
		this.localizer=localizer;
		this.breadth=parameters.getScoringBreadthType();
		this.localizationQueue=localizationQueue;
		this.minimumScore=-Log.log10(parameters.getPercolatorThreshold());
		if (parameters.getLocalizingModification().isPresent()) {
			this.localizingModification=parameters.getLocalizingModification().get();
		} else {
			throw new EncyclopediaException("You must specify a localizing modification before running localization!");
		}
	}

	@Override
	protected Nothing process() {
		EncyclopediaScorer eScorer=(EncyclopediaScorer)scorer;
		int movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		for (LibraryEntry seedEntry : super.entries) {
			ArrayList<Spectrum> stripeList=null;
			if (breadth==ScoringBreadthType.ENTIRE_RT_WINDOW) {
				stripeList=PhosphoLocalizer.getScanSubsetFromStripes(-Float.MAX_VALUE, Float.MAX_VALUE, super.stripes);
				
			} else if (breadth==ScoringBreadthType.UNCALIBRATED_20_PERCENT) {
				float duration=(super.stripes.get(super.stripes.size()-1).getScanStartTime()-super.stripes.get(0).getScanStartTime())/20.0f;
				//System.out.println(seedEntry.getScanStartTime()+" +/- "+duration); //FIXME
				stripeList=PhosphoLocalizer.getScanSubsetFromStripes(seedEntry.getScanStartTime()-duration, seedEntry.getScanStartTime()+duration, super.stripes);
				
			} else if (breadth==ScoringBreadthType.UNCALIBRATED_PEAK_WIDTH) {
				float duration=parameters.getExpectedPeakWidth();
				stripeList=PhosphoLocalizer.getScanSubsetFromStripes(seedEntry.getScanStartTime()-duration, seedEntry.getScanStartTime()+duration, super.stripes);
				
			} else if (breadth==ScoringBreadthType.RECALIBRATED_20_PERCENT||breadth==ScoringBreadthType.RECALIBRATED_PEAK_WIDTH) {
				FragmentationModel model=new FragmentationModel(seedEntry.getPeptideModSeq(), parameters.getAAConstants());
				FragmentIon[] allIons=model.getPrimaryIonObjects(parameters.getFragType(), seedEntry.getPrecursorCharge());
				float[] primary=new float[super.stripes.size()];
				for (int i=0; i<super.stripes.size(); i++) {
					Spectrum stripe=super.stripes.get(i);
					primary[i]=eScorer.score(seedEntry, stripe, allIons);
				}
				
				float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);

				float bestScore=-Float.MAX_VALUE;
				Spectrum bestStripe=null;
				for (int i=0; i<averagePrimary.length; i++) {
					if (bestScore<averagePrimary[i]) {
						bestScore=averagePrimary[i];
						bestStripe=super.stripes.get(i);
					}
				}
				
				if (breadth==ScoringBreadthType.RECALIBRATED_20_PERCENT) {
					float duration=(super.stripes.get(super.stripes.size()-1).getScanStartTime()-super.stripes.get(0).getScanStartTime())/20.0f;
					stripeList=PhosphoLocalizer.getScanSubsetFromStripes(bestStripe.getScanStartTime()-duration, bestStripe.getScanStartTime()+duration, super.stripes);
					
				} else if (breadth==ScoringBreadthType.RECALIBRATED_PEAK_WIDTH) {
					float duration=parameters.getExpectedPeakWidth();
					stripeList=PhosphoLocalizer.getScanSubsetFromStripes(bestStripe.getScanStartTime()-duration, bestStripe.getScanStartTime()+duration, super.stripes);
				}
			}

			if (stripeList==null) {
				throw new EncyclopediaException("Unexpected CASiL Scoring Breadth: "+breadth);
			}
			
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(seedEntry.getPeptideModSeq(), parameters.getAAConstants());
			
			ArrayList<String> peptideModSeqs=PhosphoPermuter.getPermutations(seedEntry.getPeptideModSeq(), localizingModification, parameters.getAAConstants());

			HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
			for (String peptideModSeq : peptideModSeqs) {
				FragmentationModel model=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
				entryMap.put(peptideModSeq, model);
			}

			ArrayList<AmbiguousPeptideModSeq> previouslyIdentified=new ArrayList<AmbiguousPeptideModSeq>();
			TIntHashSet takenRetentionTimes=new TIntHashSet();
			FragmentIonBlacklist takenIdentifiedIons=new FragmentIonBlacklist(parameters.getFragmentTolerance());
			int leftIndex=0; boolean keepLeft=true;
			int rightIndex=peptideModSeqs.size()-1; boolean keepRight=true;
			while (leftIndex<peptideModSeqs.size()&&rightIndex>=0) {
				//System.out.println("considering "+peptideModSeqs.size()+" peptideModSeqs");
				
				ArrayList<Pair<AmbiguousPeptideModSeq, FragmentIon[]>> batch=new ArrayList<Pair<AmbiguousPeptideModSeq,FragmentIon[]>>();
				
				Pair<AmbiguousPeptideModSeq, FragmentIon[]> leftPeptide;
				if (keepLeft) {
					AmbiguousPeptideModSeq leftAmbiguity=AmbiguousPeptideModSeq.getLeftAmbiguity(peptideModSeqs.get(leftIndex), localizingModification, parameters.getAAConstants());
	
					HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
					// shrink the number of unique ions subtractors to the pool of remaining sequences to the right
					for (int j=peptideModSeqs.size()-1; j>=leftIndex; j--) {
						String seq=peptideModSeqs.get(j);
						modelBatch.put(seq, entryMap.get(seq));
					}
					FragmentIon[] leftTargets=PhosphoLocalizer.getUniqueFragmentIons(leftAmbiguity.getPeptideModSeq(), seedEntry.getPrecursorCharge(), modelBatch, parameters);
					leftPeptide=new Pair<AmbiguousPeptideModSeq, FragmentIon[]>(leftAmbiguity, leftTargets);
					batch.add(leftPeptide);
					//System.out.println("ADDING LEFT: "+leftAmbiguity.getPeptideAnnotation()+" "+leftIndex+" - "+rightIndex); //FIXME
					leftIndex++;
				} else {
					leftPeptide=null;
				}
				
				Pair<AmbiguousPeptideModSeq, FragmentIon[]> rightPeptide;
				if (keepRight) {
					AmbiguousPeptideModSeq rightAmbiguity=AmbiguousPeptideModSeq.getRightAmbiguity(peptideModSeqs.get(rightIndex), localizingModification, parameters.getAAConstants());
					HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
					// shrink the number of unique ions subtractors to the pool of remaining sequences to the left
					for (int j=0; j<=rightIndex; j++) {
						String seq=peptideModSeqs.get(j);
						modelBatch.put(seq, entryMap.get(seq));
					}
					FragmentIon[] rightTargets=PhosphoLocalizer.getUniqueFragmentIons(rightAmbiguity.getPeptideModSeq(), seedEntry.getPrecursorCharge(), modelBatch, parameters);
					rightPeptide=new Pair<AmbiguousPeptideModSeq, FragmentIon[]>(rightAmbiguity, rightTargets);
					batch.add(rightPeptide);
					//System.out.println("ADDING RIGHT: "+rightAmbiguity.getPeptideAnnotation()+" "+leftIndex+" - "+rightIndex); //FIXME
					rightIndex--;
				} else {
					rightPeptide=null;
				}
				
				TIntHashSet localTakenRetentionTimes=new TIntHashSet();		
				FragmentIonBlacklist locallyTakenIdentifiedIons=new FragmentIonBlacklist(parameters.getFragmentTolerance());		
				for (Pair<AmbiguousPeptideModSeq, FragmentIon[]> pair : batch) {
					AmbiguousPeptideModSeq peptideModSeq=pair.x;
					FragmentIon[] targetIons=pair.y;
					
					//System.out.println(peptideModSeq.getPeptideAnnotation()+" "+leftIndex+" - "+rightIndex+" ("+(pair==leftPeptide?"LEFT":"RIGHT")+")"); //FIXME

					// it's ok that we don't update the targetIons based on previously IDed, for example:
					// if we know we've seen: (S[+80])SSSSK
					// and we're considering: (S[+80]S)SSSK
					// then it's ok if we use matching ions from to (S[+80])SSSSK to identify S(S[+80])SSSK
					
					// fix ambiguity based on previously identified peptides
					Optional<AmbiguousPeptideModSeq> ambiguityRemoved=peptideModSeq.removeAmbiguity(previouslyIdentified);
					if (!ambiguityRemoved.isPresent()) {
						//System.out.println("Removed ambiguity in "+peptideModSeq.getPeptideAnnotation()); //FIXME
						continue;
					}
					peptideModSeq=ambiguityRemoved.get();
					//System.out.println("considering: "+peptideModSeq.getPeptideAnnotation()); // FIXME
									
					Pair<FragmentationModel, LibraryEntry> localizedForm=seedEntry.getEntryFromNewSequence(peptideModSeq.getPeptideModSeq(), seedEntry.getAccessions(), seedEntry.isDecoy(), parameters);
					FragmentationModel localizedModel=localizedForm.x;
					LibraryEntry localizedEntry=localizedForm.y;

					AuxillaryPSMScorer auxScorer=eScorer.getAuxScorer().getEntryOptimizedScorer(localizedEntry);
					FragmentIon[] allIons=localizedModel.getPrimaryIonObjects(parameters.getFragType(), localizedEntry.getPrecursorCharge());
					
					for (int i=0; i<Math.min(localizedEntry.getMassArray().length, allIons.length); i++) {
						//System.out.println(i+") "+localizedEntry.getMassArray()[i]+"\t"+allIons[i].mass+"\t"+allIons[i].toString());
					}

					float[] primary=new float[stripeList.size()];
					for (int i=0; i<stripeList.size(); i++) {
						Spectrum stripe=stripeList.get(i);
						primary[i]=eScorer.score(localizedEntry, stripe, allIons);
					}
					
					
					float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);

					TFloatFloatHashMap scoreByRTMap=new TFloatFloatHashMap();
					ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
					for (int i=0; i<averagePrimary.length; i++) {
						
						goodStripes.add(new ScoredIndex(averagePrimary[i], i));
						scoreByRTMap.put(stripeList.get(i).getScanStartTime(), primary[i]);
					}
					Collections.sort(goodStripes);

					EValueCalculator calculator=new EValueCalculator(scoreByRTMap);

					PeptideScoringResult result=new PeptideScoringResult(localizedEntry);
					// find the best stripe and localize near it
					int identifiedPeaks=0;
					for (int i=goodStripes.size()-1; i>=0; i--) {
						int index=goodStripes.get(i).y;
						Spectrum stripe=stripeList.get(index);
						int stripeRTIndex=getStripeRTIndex(stripe);
						if (takenRetentionTimes.contains(stripeRTIndex)) {
							continue;
							
						} else {
							ArrayList<Spectrum> stripeSubset=PhosphoLocalizer.getScanSubsetFromStripes(stripe.getScanStartTime()-parameters.getExpectedPeakWidth(), stripe.getScanStartTime()+parameters.getExpectedPeakWidth(), stripes);
							Triplet<ModificationLocalizationData, Stripe, Range> locData=calculateLocalizationScoring(minimumScore, parameters, dutyCycle, localizer, seedEntry, peptideModSeq, targetIons, allIons, takenIdentifiedIons, stripeSubset);

							ModificationLocalizationData data=locData.x;
							Stripe apex=locData.y;
							Range peakRange=locData.z;
							
							float[] auxScoreArray=auxScorer.score(localizedEntry, apex, predictedIsotopeDistribution, precursors);

							float score=eScorer.score(localizedEntry, apex, allIons);
							float evalue=calculator.getNegLog10EValue(score);
							if (Float.isNaN(evalue)) {
								evalue=-1.0f;
							}

							result.addStripe(score, General.concatenate(auxScoreArray, evalue, data.getLocalizationScore()), apex);
							resultsQueue.add(result);
							
							//System.out.println("\tlocalization:"+data.isSiteSpecific()+"\t"+data.getLocalizationPeptideModSeq().getPeptideAnnotation()+"\t"+apex.getScanStartTime()+"\t"+data.getLocalizationScore()+"\t"+FragmentIon.toArchiveString(data.getLocalizingIons()));
							if (!localizedEntry.isDecoy()) {
								// don't bother logging decoys
								localizationQueue.add(data);
							}
							
							if (!data.isSiteSpecific()) {
								if (pair==leftPeptide) {
									keepLeft=false;
								} else {
									keepRight=false;
								}
							} else {								
								// allows searching beyond this mod (only if we localize it to a RT)
								previouslyIdentified.add(peptideModSeq);
								entryMap.remove(peptideModSeq); 
							
								// block +/- a peakWidth window
								int removalIndex=index;
								while (removalIndex>=0&&peakRange.contains(stripeList.get(removalIndex).getScanStartTime())) {
									localTakenRetentionTimes.add(removalIndex);
									removalIndex--;
								}
								removalIndex=index+1;
								while (removalIndex<stripeList.size()&&peakRange.contains(stripeList.get(removalIndex).getScanStartTime())) {
									localTakenRetentionTimes.add(removalIndex);
									removalIndex++;
								}
	
								for (FragmentIon target : allIons) {
									locallyTakenIdentifiedIons.addIonToBlacklist(target.mass, peakRange);
								}
							}
							
							if (identifiedPeaks>peaksKept) {
								// keep N+1 peaks
								break;
							}

							identifiedPeaks++;
						}
					}
				}
				if (!keepLeft&&!keepRight) {
					// wasn't able to localize either side, so can't keep going
					//System.out.println("BREAKING BATCH!"); //FIXME
					break;
					
				} else {
					// Note, it's ok if opposite sides (in the pairs) go head to head at the same RT, since they implicitly consider different fragment ions
					takenRetentionTimes.addAll(localTakenRetentionTimes);
					takenIdentifiedIons.addIonsToBlacklist(locallyTakenIdentifiedIons);
				}
			}
		}
		return Nothing.NOTHING;
	}

	/**
	 * straight truncate, NOTE, this means we don't have peak resolution less than 1 second
	 * @param stripe
	 * @return
	 */
	private int getStripeRTIndex(Spectrum stripe) {
		return (int)stripe.getScanStartTime();
	}

	public static Triplet<ModificationLocalizationData, Stripe, Range> calculateLocalizationScoring(float minimumScore, SearchParameters parameters, float dutyCycle, PhosphoLocalizer localizer, LibraryEntry seedEntry, AmbiguousPeptideModSeq peptideModSeq, FragmentIon[] targetIons, FragmentIon[] allIons, FragmentIonBlacklist takenIdentifiedIons, ArrayList<Spectrum> stripeSubset) {
		int targetNumFragments=Math.max(parameters.getMinNumOfQuantitativePeaks(), 3);
									
		double[] targetIonsMasses=FragmentIon.getMasses(targetIons);
		float[] frequencies=localizer.getBackground().getFrequencies(targetIonsMasses, seedEntry.getPrecursorMZ(), parameters.getFragmentTolerance());
		
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" --> "+General.toString(targetIonsMasses));
		
		float[] negLogProbsSiteSpecific=new float[stripeSubset.size()];
		for (int k=0; k<stripeSubset.size(); k++) {
			negLogProbsSiteSpecific[k]=PhosphoLocalizer.score(parameters, targetIonsMasses, targetIons, frequencies, stripeSubset.get(k), true);
		}
		float bestLocalizationScore=-Float.MAX_VALUE;
		Stripe apex=null;
		negLogProbsSiteSpecific=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, Math.round(parameters.getExpectedPeakWidth()/dutyCycle));//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);
		for (int k=0; k<stripeSubset.size(); k++) {
			if (bestLocalizationScore<negLogProbsSiteSpecific[k]) {
				bestLocalizationScore=negLogProbsSiteSpecific[k];
				apex=(Stripe)stripeSubset.get(k);
			}
		}
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" --> ("+targetIons.length+") "+bestLocalizationScore+" localization score at "+(apex.getScanStartTime()/60f)+" min"); //FIXME 

		boolean wasLocalized=false;
		float localizationIntensity=0.0f;
		float totalIntensity=0.0f;
		int numIdentificationPeaks=0;
		int numberOfMods=PeptideUtils.getNumberOfMods(peptideModSeq.getPeptideModSeq(), localizer.getModification().getNominalMass());
		ArrayList<FragmentIon> wellShapedIons=new ArrayList<FragmentIon>();
		
		Range peakRange=new Range(apex.getScanStartTime(), apex.getScanStartTime());
		if (bestLocalizationScore>=minimumScore) {
			// generate quant data from localizing ions only
			TransitionRefinementData quantData=localizer.quantifyPeptide(peptideModSeq.getPeptideModSeq(), seedEntry.getPrecursorCharge(), targetIons, apex.getScanStartTime(),
					stripeSubset, takenIdentifiedIons, Optional.ofNullable((float[]) null));
			if (quantData!=null) {
				peakRange=quantData.getRange();
				float[] intensities=quantData.getIntegrationArray();
				float[] correlations=quantData.getCorrelationArray();
				FragmentIon[] consideredIons=quantData.getFragmentMassArray();
				for (int i=0; i<consideredIons.length; i++) {
					wellShapedIons.add(consideredIons[i]);
					localizationIntensity+=intensities[i];
				}

				// calculate quant data for all ions
				float[] medianChromatogram=quantData.getMedianChromatogram();
				TransitionRefinementData allQuantData=localizer.quantifyPeptide(peptideModSeq.getPeptideModSeq(), seedEntry.getPrecursorCharge(), allIons, apex.getScanStartTime(),
						stripeSubset, takenIdentifiedIons, Optional.of(medianChromatogram));
				if (allQuantData!=null) {
					peakRange=allQuantData.getRange();
					intensities=allQuantData.getIntegrationArray();
					correlations=allQuantData.getCorrelationArray();
					for (int i=0; i<correlations.length; i++) {
						if (correlations[i]>=TransitionRefiner.identificationCorrelationThreshold) {
							numIdentificationPeaks++;
							totalIntensity+=intensities[i];
						}
					}
					wasLocalized=wellShapedIons.size()>0&&numIdentificationPeaks>=targetNumFragments&&AmbiguousPeptideModSeq.isLocalized(peptideModSeq, localizer.getModification());
					//System.out.println("\tLocalized "+wasLocalized+" for "+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score, "+numIdentificationPeaks+"/"+correlations.length+" peaks)"); // FIXME
				}
			}
		}
		
		ModificationLocalizationData modData=new ModificationLocalizationData(peptideModSeq, apex.getScanStartTime(), bestLocalizationScore, numberOfMods, wasLocalized, wellShapedIons.toArray(new FragmentIon[wellShapedIons.size()]), localizationIntensity, totalIntensity);
		return new Triplet<ModificationLocalizationData, Stripe, Range>(modData, apex, peakRange);
	}
	
}
