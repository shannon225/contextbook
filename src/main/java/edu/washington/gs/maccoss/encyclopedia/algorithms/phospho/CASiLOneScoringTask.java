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
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
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
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.list.array.TDoubleArrayList;
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
				FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, parameters.getAAConstants());
				entryMap.put(peptideModSeq, model);
			}

			// generate an ion map for all forms
			HashMap<String, double[]> ionsByPeptide=new HashMap<>();
			for (Entry<String, FragmentationModel> modelEntry : entryMap.entrySet()) {
				ionsByPeptide.put(modelEntry.getKey(), modelEntry.getValue().getPrimaryIons(parameters.getFragType(), firstEntry.getPrecursorCharge(), true));
			}

			ArrayList<AmbiguousPeptideModSeq> previouslyIdentified=new ArrayList<AmbiguousPeptideModSeq>();
			TIntHashSet takenRetentionTimes=new TIntHashSet();
			FragmentIonBlacklist takenIdentifiedIons=new FragmentIonBlacklist(parameters.getFragmentTolerance());
			int leftIndex=0; boolean keepLeft=true;
			int rightIndex=peptideModSeqs.size()-1; boolean keepRight=true;
			while ((keepLeft&&leftIndex<peptideModSeqs.size())||(keepRight&&rightIndex>=0)) {
				//System.out.println("considering "+peptideModSeqs.size()+" peptideModSeqs");
				
				ArrayList<Pair<AmbiguousPeptideModSeq, FragmentIon[]>> batch=new ArrayList<Pair<AmbiguousPeptideModSeq,FragmentIon[]>>();
				
				Pair<AmbiguousPeptideModSeq, FragmentIon[]> leftPeptide = getLeftPeptide(precursorCharge, peptideModSeqs, entryMap, leftIndex, keepLeft);
				if (leftPeptide!=null) {
					batch.add(leftPeptide);
					//System.out.println("ADDING LEFT: "+leftAmbiguity.getPeptideAnnotation()+" "+leftIndex+" - "+rightIndex); //FIXME
					leftIndex++;
				}
				
				Pair<AmbiguousPeptideModSeq, FragmentIon[]> rightPeptide = getRightPeptide(precursorCharge, peptideModSeqs, entryMap, rightIndex, keepRight);
				if (rightPeptide!=null) {
					batch.add(rightPeptide);
					//System.out.println("ADDING RIGHT: "+rightAmbiguity.getPeptideAnnotation()+" "+leftIndex+" - "+rightIndex); //FIXME
					rightIndex--;
				}
				
				TIntHashSet localTakenRetentionTimes=new TIntHashSet();		
				FragmentIonBlacklist locallyTakenIdentifiedIons=new FragmentIonBlacklist(parameters.getFragmentTolerance());		
				for (Pair<AmbiguousPeptideModSeq, FragmentIon[]> pair : batch) {
					boolean localizeLeftSide = pair==leftPeptide;
					
					AmbiguousPeptideModSeq peptideModSeq=pair.x;
					FragmentIon[] targetIons=pair.y;
					
					//System.out.println(peptideModSeq.getPeptideAnnotation()+" "+leftIndex+" - "+rightIndex+" ("+(pair==leftPeptide?"LEFT":"RIGHT")+")"); //FIXME

					// it's ok that we don't update the targetIons based on previously IDed, for example:
					// if we know we've seen: (S[+80])SSSSK
					// and we're considering: (S[+80]S)SSSK
					// then it's ok if we use matching ions from to (S[+80])SSSSK to identify S(S[+80])SSSK
					
					// fix ambiguity based on previously identified peptides
					Optional<AmbiguousPeptideModSeq> ambiguityRemoved=peptideModSeq.removeAmbiguity(localizingModification, previouslyIdentified);
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
					FragmentIon[] allIons=localizedModel.getPrimaryIonObjects(parameters.getFragType(), localizedEntry.getPrecursorCharge(), false, true);
					
					float[] primary=new float[stripeList.size()];
					for (int i=0; i<stripeList.size(); i++) {
						Spectrum stripe=stripeList.get(i);

						int stripeRTIndex=getStripeRTIndex(stripe);
						if (takenRetentionTimes.contains(stripeRTIndex)) {
							primary[i]=-1.0f;
						} else {
							primary[i]=eScorer.score(localizedEntry, stripe, allIons);
						}
					}
					
					float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);
					averagePrimary=General.subtract(averagePrimary, QuickMedian.median(averagePrimary.clone()));
					
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
							boolean wasPreviouslyObservedSite=localTakenRetentionTimes.contains(stripeRTIndex); // implicitly: ||takenRetentionTimes.contains(stripeRTIndex)
							ArrayList<Spectrum> stripeSubset=PhosphoLocalizer.getScanSubsetFromStripes(stripe.getScanStartTime()-parameters.getExpectedPeakWidth(), stripe.getScanStartTime()+parameters.getExpectedPeakWidth(), stripes);
							Triplet<ModificationLocalizationData, Stripe, Range> locData=calculateLocalizationScoring(wasPreviouslyObservedSite, minimumScore, parameters, dutyCycle, localizer, localizedEntry, peptideModSeq, targetIons, allIons, takenIdentifiedIons, stripeSubset);
							while (!locData.x.isLocalized()) {
								if (localizeLeftSide) {
									leftPeptide=getLeftPeptide(precursorCharge, peptideModSeqs, entryMap, leftIndex, keepLeft);
									if (leftPeptide==null) break;
									//System.out.println("Couldn't localize, expanding further left ("+leftIndex+")...");
									leftIndex++;
									peptideModSeq=leftPeptide.x;
									targetIons=leftPeptide.y;
									locData=calculateLocalizationScoring(wasPreviouslyObservedSite, minimumScore, parameters, dutyCycle, localizer, localizedEntry, peptideModSeq, targetIons, allIons, takenIdentifiedIons, stripeSubset);
								} else {
									rightPeptide=getRightPeptide(precursorCharge, peptideModSeqs, entryMap, rightIndex, keepRight);
									if (rightPeptide==null) break;
									//System.out.println("Couldn't localize, expanding further right ("+rightIndex+")...");
									rightIndex--;
									peptideModSeq=rightPeptide.x;
									targetIons=rightPeptide.y;
									locData=calculateLocalizationScoring(wasPreviouslyObservedSite, minimumScore, parameters, dutyCycle, localizer, localizedEntry, peptideModSeq, targetIons, allIons, takenIdentifiedIons, stripeSubset);
								}
							}
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
							
							//System.out.println("\tlocalization:"+data.isSiteSpecific()+"/"+data.isLocalized()+"/"+(data.getLocalizationPeptideModSeq().getNumModifiableSites()+"=="+data.getNumberOfMods())+"\t"+data.getLocalizationPeptideModSeq().getPeptideAnnotation()+"\t"+apex.getScanStartTime()+"\t"+data.getLocalizationScore()+"\t"+FragmentIon.toArchiveString(data.getLocalizingIons()));
							if (!localizedEntry.isDecoy()) {
								// don't bother logging decoys
								localizationQueue.add(data);
							}
							
							if (!data.isSiteSpecific()) {
								if (localizeLeftSide) {
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
									localTakenRetentionTimes.add(getStripeRTIndex(stripeList.get(removalIndex)));
									removalIndex--;
								}
								removalIndex=index+1;
								while (removalIndex<stripeList.size()&&peakRange.contains(stripeList.get(removalIndex).getScanStartTime())) {
									localTakenRetentionTimes.add(getStripeRTIndex(stripeList.get(removalIndex)));
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

	private Pair<AmbiguousPeptideModSeq, FragmentIon[]> getRightPeptide(byte precursorCharge,
			ArrayList<String> peptideModSeqs, HashMap<String, FragmentationModel> entryMap, int rightIndex,
			boolean keepRight) {
		Pair<AmbiguousPeptideModSeq, FragmentIon[]> rightPeptide;
		if (keepRight&&rightIndex>=0) {
			AmbiguousPeptideModSeq rightAmbiguity=AmbiguousPeptideModSeq.getRightAmbiguity(peptideModSeqs.get(rightIndex), localizingModification, parameters.getAAConstants(), "");
			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of remaining sequences to the left
			for (int j=0; j<=rightIndex; j++) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			FragmentIon[] rightTargets=PhosphoLocalizer.getUniqueFragmentIons(rightAmbiguity.getPeptideModSeq(), precursorCharge, modelBatch, parameters);
			rightPeptide=new Pair<AmbiguousPeptideModSeq, FragmentIon[]>(rightAmbiguity, rightTargets);
			
		} else {
			rightPeptide=null;
		}
		return rightPeptide;
	}

	private Pair<AmbiguousPeptideModSeq, FragmentIon[]> getLeftPeptide(byte precursorCharge,
			ArrayList<String> peptideModSeqs, HashMap<String, FragmentationModel> entryMap, int leftIndex,
			boolean keepLeft) {
		Pair<AmbiguousPeptideModSeq, FragmentIon[]> leftPeptide;
		if (keepLeft&&leftIndex<peptideModSeqs.size()) {
			AmbiguousPeptideModSeq leftAmbiguity=AmbiguousPeptideModSeq.getLeftAmbiguity(peptideModSeqs.get(leftIndex), localizingModification, parameters.getAAConstants(), "");

			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of remaining sequences to the right
			for (int j=peptideModSeqs.size()-1; j>=leftIndex; j--) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			FragmentIon[] leftTargets=PhosphoLocalizer.getUniqueFragmentIons(leftAmbiguity.getPeptideModSeq(), precursorCharge, modelBatch, parameters);
			leftPeptide=new Pair<AmbiguousPeptideModSeq, FragmentIon[]>(leftAmbiguity, leftTargets);
			
		} else {
			leftPeptide=null;
		}
		return leftPeptide;
	}
	
	static class LocalizedForm {
		private final FragmentationModel localizedModel;
		private final LibraryEntry localizedEntry;
		private final ArrayList<Spectrum> scansToConsider;
		public LocalizedForm(FragmentationModel localizedModel, LibraryEntry localizedEntry, ArrayList<Spectrum> scansToConsider) {
			this.localizedModel=localizedModel;
			this.localizedEntry=localizedEntry;
			this.scansToConsider=scansToConsider;
		}
	}
	
	public static LocalizedForm getLocalizedForm(AmbiguousPeptideModSeq targetPeptide, byte charge, HashMap<String, FragmentationModel> modelMap, HashMap<String, double[]> ionsByPeptide, HashMap<LibraryEntry, ArrayList<Spectrum>> scansByEntry, SearchParameters parameters) {
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
			// This can happen when we're localizing around an n-term acetyl that's been rearranged incorrectly by reversing
			//Logger.errorLine("Missing target ions for: "+targetPeptideModSeq+", Found ions for: "+General.toString(ionsByPeptide.keySet())+", skipping form");
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
				//Logger.errorLine("Missing real ions for: "+realEntry.getAccuratePeptideModSeq(parameters.getAAConstants())+", Found ions for: "+General.toString(ionsByPeptide.keySet()));
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
	
	private static int getNumberOfMatchingIons(double[] a, double[] b, MassTolerance tolerance) {
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
			FragmentationModel model=PeptideUtils.getPeptideModel(seedEntry.getPeptideModSeq(), parameters.getAAConstants());
			FragmentIon[] allIons=model.getPrimaryIonObjects(parameters.getFragType(), seedEntry.getPrecursorCharge(), false, true);
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
	 * straight truncate, NOTE, this means we don't have peak resolution less than 0.1 second
	 * @param stripe
	 * @return
	 */
	private int getStripeRTIndex(Spectrum stripe) {
		return (int)(stripe.getScanStartTime()*10);
	}

	public static Triplet<ModificationLocalizationData, Stripe, Range> calculateLocalizationScoring(boolean wasPreviouslyObservedSite, float minimumScore, SearchParameters parameters, float dutyCycle, PhosphoLocalizer localizer, LibraryEntry localizedEntry, AmbiguousPeptideModSeq peptideModSeq, FragmentIon[] targetIons, FragmentIon[] allIons, FragmentIonBlacklist takenIdentifiedIons, ArrayList<Spectrum> stripeSubset) {
		Pair<Stripe, Float> bestLocalizedStripe = getBestLocalizationStripe(parameters, dutyCycle, localizer, localizedEntry, targetIons, stripeSubset);
		return generateLocalizationData(wasPreviouslyObservedSite, minimumScore, parameters, localizer, localizedEntry,
				peptideModSeq, targetIons, allIons, takenIdentifiedIons, stripeSubset, bestLocalizedStripe);
	}

	public static Triplet<ModificationLocalizationData, Stripe, Range> generateLocalizationData(
			boolean wasPreviouslyObservedSite, float minimumScore, SearchParameters parameters,
			PhosphoLocalizer localizer, LibraryEntry localizedEntry, AmbiguousPeptideModSeq peptideModSeq,
			FragmentIon[] targetIons, FragmentIon[] allIons, FragmentIonBlacklist takenIdentifiedIons,
			ArrayList<Spectrum> stripeSubset, Pair<Stripe, Float> bestLocalizedStripe) {
		float bestLocalizationScore=bestLocalizedStripe.y;
		Stripe apex=bestLocalizedStripe.x;
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" --> ("+targetIons.length+") "+bestLocalizationScore+" localization score at "+(apex.getScanStartTime()/60f)+" min"); //FIXME 
		
		int targetNumFragments=Math.max(parameters.getMinNumOfQuantitativePeaks(), 3);

		boolean isLocalized=false;
		boolean isSiteSpecific=false;
		boolean isCompletelyAmbiguous=false;
		float localizationIntensity=0.0f;
		float totalIntensity=0.0f;
		float numIdentificationPeaks=0.0f;
		int numConsideredPeaks=0;
		int numberOfMods=PeptideUtils.getNumberOfMods(peptideModSeq.getPeptideModSeq(), localizer.getModification().getNominalMass());
		ArrayList<FragmentIon> wellShapedIons=new ArrayList<FragmentIon>();

        Range peakRange=new Range(apex.getScanStartTime(), apex.getScanStartTime());
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score)\tNOT GOOD ENOUGH");
		if (bestLocalizationScore>=minimumScore||(bestLocalizationScore>0&&!wasPreviouslyObservedSite)) {
			// generate quant data from localizing ions only
			TransitionRefinementData quantData=localizer.quantifyPeptide(parameters, peptideModSeq.getPeptideModSeq(), localizedEntry.getPrecursorCharge(), targetIons, apex.getScanStartTime(),
					stripeSubset, takenIdentifiedIons, Optional.ofNullable((float[]) null));
			//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score)\ta)"+(quantData!=null));
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
				TransitionRefinementData allQuantData=localizer.quantifyPeptide(parameters, peptideModSeq.getPeptideModSeq(), localizedEntry.getPrecursorCharge(), allIons, apex.getScanStartTime(),
						stripeSubset, takenIdentifiedIons, Optional.of(medianChromatogram));
				//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score)\tb)"+(allQuantData!=null));
				if (allQuantData!=null&&allQuantData.getMassArray().isPresent()) {
					peakRange=allQuantData.getRange();
					intensities=allQuantData.getIntegrationArray();
					correlations=allQuantData.getCorrelationArray();
					
					TDoubleArrayList matchingPeaks=new TDoubleArrayList();
					for (int i=0; i<correlations.length; i++) {
						numConsideredPeaks++;
						/*if (correlations[i]>=TransitionRefiner.quantitativeCorrelationThreshold&&intensities[i]>0.0f) {
							numIdentificationPeaks+=1.0f;
							totalIntensity+=intensities[i];
						} else if (correlations[i]>=TransitionRefiner.identificationCorrelationThreshold&&intensities[i]>0.0f) {
							numIdentificationPeaks+=0.1f;
						} else if (intensities[i]>0.0f) {
							numIdentificationPeaks+=0.1f;
						}*/
						if (correlations[i]>0.0f&&intensities[i]>0.0f) {
							numIdentificationPeaks+=correlations[i]*correlations[i];
						}
						if (intensities[i]>0.0f) {
							totalIntensity+=intensities[i];
							matchingPeaks.add(allIons[i].mass);
						}
						//System.out.println(localizedEntry.getPeptideModSeq()+", "+numIdentificationPeaks+") "+targetIonSet.contains(allIons[i])+", "+allIons[i].toString()+", "+allIons[i].mass+"\t"+intensities[i]+"\t"+correlations[i]);
					}
					// recalculate localization scores using only the peaks are remotely ok
					// FIXME
//					double[] targetIonsMasses=FragmentIon.getMasses(targetIons);
//					float[] frequencies=localizer.getBackground().getFrequencies(targetIonsMasses, localizedEntry.getPrecursorMZ(), parameters.getFragmentTolerance());
//					bestLocalizationScore=PhosphoLocalizer.score(parameters, targetIonsMasses, targetIons, frequencies, matchingPeaks.toArray());

					//System.out.println(localizedEntry.getPeptideModSeq()+" Final: "+bestLocalizationScore);
					
					isCompletelyAmbiguous=AmbiguousPeptideModSeq.isCompletelyAmbiguous(peptideModSeq, localizer.getModification());
					isLocalized=bestLocalizationScore>=minimumScore&&wellShapedIons.size()>0&&numIdentificationPeaks>targetNumFragments&&!isCompletelyAmbiguous;
					isSiteSpecific=isLocalized&&AmbiguousPeptideModSeq.isSiteSpecific(peptideModSeq, localizer.getModification());
					//System.out.println("\tLocalized "+isLocalized+"/"+isSiteSpecific+" for "+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score vs "+minimumScore+"minimum, "+numIdentificationPeaks+"/"+correlations.length+" peaks)"); // FIXME
				}
			}
		}
		//System.out.println("CASiL: "+peptideModSeq.getPeptideModSeq()+" ("+targetIons.length+"/"+allIons.length+"/"+numConsideredPeaks+")\t"+bestLocalizationScore+","+isCompletelyAmbiguous+","+isLocalized+","+isSiteSpecific+","+numIdentificationPeaks);
		
		ModificationLocalizationData modData=new ModificationLocalizationData(peptideModSeq, apex.getScanStartTime(), bestLocalizationScore, numberOfMods, isSiteSpecific, isLocalized, isCompletelyAmbiguous, wellShapedIons.toArray(new FragmentIon[wellShapedIons.size()]), localizationIntensity, totalIntensity);
		return new Triplet<ModificationLocalizationData, Stripe, Range>(modData, apex, peakRange);
	}

	public static Pair<Stripe, Float> getBestLocalizationStripe(SearchParameters parameters, float dutyCycle,
			PhosphoLocalizer localizer, LibraryEntry localizedEntry, FragmentIon[] targetIons,
			ArrayList<Spectrum> stripeSubset) {
		double[] targetIonsMasses=FragmentIon.getMasses(targetIons);
		float[] frequencies=localizer.getBackground().getFrequencies(targetIonsMasses, localizedEntry.getPrecursorMZ(), parameters.getFragmentTolerance());
		
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" --> "+General.toString(targetIonsMasses));

		float[] negLogProbsSiteSpecific=new float[stripeSubset.size()];
		for (int k=0; k<stripeSubset.size(); k++) {
			negLogProbsSiteSpecific[k]=PhosphoLocalizer.score(parameters, targetIonsMasses, targetIons, frequencies, stripeSubset.get(k), true);
		}
		float bestScore=-Float.MAX_VALUE;
		Stripe apexStripe=null;
		negLogProbsSiteSpecific=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, Math.round(parameters.getExpectedPeakWidth()/dutyCycle));//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);
		for (int k=0; k<stripeSubset.size(); k++) {
			if (bestScore<negLogProbsSiteSpecific[k]) {
				bestScore=negLogProbsSiteSpecific[k];
				apexStripe=(Stripe)stripeSubset.get(k);
			}
		}
		
		Pair<Stripe, Float> bestLocalizedStripe=new Pair<Stripe, Float>(apexStripe, bestScore);
		return bestLocalizedStripe;
	}
	
}
