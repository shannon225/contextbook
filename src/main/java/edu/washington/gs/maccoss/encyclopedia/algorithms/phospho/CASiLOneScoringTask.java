package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
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
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
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
		HashMap<String, ArrayList<LibraryEntry>> entriesBySequence=new HashMap<>();
		for (LibraryEntry entry : super.entries) {
			String seq=entry.getPeptideSeq();
			ArrayList<LibraryEntry> list=entriesBySequence.get(seq);
			if (list==null) {
				list=new ArrayList<>();
				entriesBySequence.put(seq, list);
			}
			list.add(entry);
		}// FIXME FLESH OUT REPEATED ENTRIES FOR LOCALIZATION
		
		for (ArrayList<LibraryEntry> seedEntries : entriesBySequence.values()) {
			LibraryEntry firstEntry=seedEntries.get(0);
			HashMap<LibraryEntry, ArrayList<Spectrum>> scansByEntry=new HashMap<>();
			for (LibraryEntry libraryEntry : seedEntries) {
				scansByEntry.put(libraryEntry, getSpectra(eScorer, movingAverageLength, libraryEntry));
			}
			
			byte precursorCharge=firstEntry.getPrecursorCharge();
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(firstEntry.getPeptideModSeq(), parameters.getAAConstants());
			ArrayList<String> peptideModSeqs=PhosphoPermuter.getPermutations(firstEntry.getPeptideModSeq(), localizingModification, parameters.getAAConstants());

			HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
			for (String peptideModSeq : peptideModSeqs) {
				FragmentationModel model=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
				entryMap.put(peptideModSeq, model);
			}

			// generate an ion map for all forms
			HashMap<String, double[]> ionsByPeptide=new HashMap<>();
			for (Entry<String, FragmentationModel> modelEntry : entryMap.entrySet()) {
				ionsByPeptide.put(modelEntry.getKey(), modelEntry.getValue().getPrimaryIons(parameters.getFragType(), firstEntry.getPrecursorCharge()));
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
					FragmentIon[] leftTargets=PhosphoLocalizer.getUniqueFragmentIons(leftAmbiguity.getPeptideModSeq(), precursorCharge, modelBatch, parameters);
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
					FragmentIon[] rightTargets=PhosphoLocalizer.getUniqueFragmentIons(rightAmbiguity.getPeptideModSeq(), precursorCharge, modelBatch, parameters);
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
									
					LocalizedForm localizedForm=getLocalizedForm(peptideModSeq, firstEntry.getPrecursorCharge(), entryMap, ionsByPeptide, scansByEntry, parameters);
					if (localizedForm==null) {
						continue;
					}
					
					FragmentationModel localizedModel=localizedForm.localizedModel;
					LibraryEntry localizedEntry=localizedForm.localizedEntry;
					ArrayList<Spectrum> stripeList=localizedForm.scansToConsider;

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
							Triplet<ModificationLocalizationData, Stripe, Range> locData=calculateLocalizationScoring(minimumScore, parameters, dutyCycle, localizer, localizedEntry, peptideModSeq, targetIons, allIons, takenIdentifiedIons, stripeSubset);

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
	
	class LocalizedForm {
		private final FragmentationModel localizedModel;
		private final LibraryEntry localizedEntry;
		private final ArrayList<Spectrum> scansToConsider;
		public LocalizedForm(FragmentationModel localizedModel, LibraryEntry localizedEntry, ArrayList<Spectrum> scansToConsider) {
			this.localizedModel=localizedModel;
			this.localizedEntry=localizedEntry;
			this.scansToConsider=scansToConsider;
		}
	}
	
	public LocalizedForm getLocalizedForm(AmbiguousPeptideModSeq targetPeptide, byte charge, HashMap<String, FragmentationModel> modelMap, HashMap<String, double[]> ionsByPeptide, HashMap<LibraryEntry, ArrayList<Spectrum>> scansByEntry, SearchParameters parameters) {
		String targetPeptideModSeq=targetPeptide.getPeptideModSeq();
		FragmentationModel targetModel=modelMap.get(targetPeptideModSeq);

		// first check if we have real data to cover this form
		for (Entry<LibraryEntry, ArrayList<Spectrum>> realDataEntry : scansByEntry.entrySet()) {
			LibraryEntry realEntry=realDataEntry.getKey();
			if (targetPeptideModSeq.equals(realEntry.getAccuratePeptideModSeq(parameters.getAAConstants()))) {
				return new LocalizedForm(targetModel, realEntry, realDataEntry.getValue());
			}
		}

		int bestNumMatching=-1;
		LibraryEntry bestRealEntry=null;
		ArrayList<Spectrum> bestScans=null;
		
		double[] targetIons=ionsByPeptide.get(targetPeptideModSeq);
		if (targetIons==null) {
			Logger.errorLine("Missing target ions for: "+targetPeptideModSeq+", Found ions for: "+General.toString(ionsByPeptide.keySet())+", skipping form");
			return null;
		}
		for (Entry<LibraryEntry, ArrayList<Spectrum>> realDataEntry : scansByEntry.entrySet()) {
			LibraryEntry realEntry=realDataEntry.getKey();
			FragmentationModel realModel=modelMap.get(realEntry.getAccuratePeptideModSeq(parameters.getAAConstants()));
			if (realModel==null) {
				// these will be associated with oxidation or other peptides within the same window (and shouldn't be used to model fragmentation)
				continue;
			}
			
			double[] realIons=ionsByPeptide.get(realEntry.getAccuratePeptideModSeq(parameters.getAAConstants()));
			if (realIons==null) {
				Logger.errorLine("Missing real ions for: "+realEntry.getAccuratePeptideModSeq(parameters.getAAConstants())+", Found ions for: "+General.toString(ionsByPeptide.keySet()));
				continue;				
			}
			int numMatching=getNumberOfMatchingIons(targetIons, realIons, parameters.getFragmentTolerance());
			
			if (numMatching>bestNumMatching) {
				bestNumMatching=numMatching;
				bestRealEntry=realEntry;
				bestScans=realDataEntry.getValue();
			}
		}

		Pair<FragmentationModel, LibraryEntry> localizedForm=bestRealEntry.getEntryFromNewSequence(targetPeptideModSeq, bestRealEntry.getAccessions(), bestRealEntry.isDecoy(), parameters);
		return new LocalizedForm(localizedForm.x, localizedForm.y, bestScans);
	}
	
	private int getNumberOfMatchingIons(double[] a, double[] b, MassTolerance tolerance) {
		if (a==null||b==null) return 0;
		double[] x,y;
		if (a.length>b.length) {
			y=a;
			x=b;
		} else {
			y=b;
			x=a;
		}
		int matches=0;
		for (int i=0; i<x.length; i++) {
			if (tolerance.getIndex(y, x[i]).isPresent()) {
				matches++;
			}
		}
		return matches;
	}

	private ArrayList<Spectrum> getSpectra(EncyclopediaScorer eScorer, int movingAverageLength, LibraryEntry seedEntry) {
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
		return stripeList;
	}

	/**
	 * straight truncate, NOTE, this means we don't have peak resolution less than 1 second
	 * @param stripe
	 * @return
	 */
	private int getStripeRTIndex(Spectrum stripe) {
		return (int)stripe.getScanStartTime();
	}

	public static Triplet<ModificationLocalizationData, Stripe, Range> calculateLocalizationScoring(float minimumScore, SearchParameters parameters, float dutyCycle, PhosphoLocalizer localizer, LibraryEntry localizedEntry, AmbiguousPeptideModSeq peptideModSeq, FragmentIon[] targetIons, FragmentIon[] allIons, FragmentIonBlacklist takenIdentifiedIons, ArrayList<Spectrum> stripeSubset) {
		int targetNumFragments=Math.max(parameters.getMinNumOfQuantitativePeaks(), 3);
									
		double[] targetIonsMasses=FragmentIon.getMasses(targetIons);
		float[] frequencies=localizer.getBackground().getFrequencies(targetIonsMasses, localizedEntry.getPrecursorMZ(), parameters.getFragmentTolerance());
		
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
			TransitionRefinementData quantData=localizer.quantifyPeptide(peptideModSeq.getPeptideModSeq(), localizedEntry.getPrecursorCharge(), targetIons, apex.getScanStartTime(),
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
				TransitionRefinementData allQuantData=localizer.quantifyPeptide(peptideModSeq.getPeptideModSeq(), localizedEntry.getPrecursorCharge(), allIons, apex.getScanStartTime(),
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
