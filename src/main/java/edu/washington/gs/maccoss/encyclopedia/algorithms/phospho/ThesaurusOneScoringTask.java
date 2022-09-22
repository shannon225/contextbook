package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Arrays;
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
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class ThesaurusOneScoringTask extends AbstractLibraryScoringTask {
	
	public static final float MAXIMUM_DELTA_FROM_BEST_SCORE = 0.9f;
	private final PhosphoLocalizer localizer;
	private final float dutyCycle;
	private final ScoringBreadthType breadth;
	private final PeptideModification localizingModification;
	private final BlockingQueue<ModificationLocalizationData> localizationQueue;
	private final float minimumScore;
	private final ThesaurusSearchParameters thesaurusParameters;
	
	public ThesaurusOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, float dutyCycle, PrecursorScanMap precursors, 
			PhosphoLocalizer localizer, BlockingQueue<AbstractScoringResult> resultsQueue, BlockingQueue<ModificationLocalizationData> localizationQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.dutyCycle=dutyCycle;
		this.localizer=localizer;
		this.breadth=parameters.getScoringBreadthType();
		this.localizationQueue=localizationQueue;
		this.minimumScore=2.0f;//-Log.log10(parameters.getPercolatorThreshold());
		if (parameters.getLocalizingModification().isPresent()) {
			this.localizingModification=parameters.getLocalizingModification().get();
		} else {
			throw new EncyclopediaException("You must specify a localizing modification before running localization!");
		}
		if (parameters instanceof ThesaurusSearchParameters) {
			this.thesaurusParameters=(ThesaurusSearchParameters)parameters;
		} else {
			throw new EncyclopediaException("You must specify a Thesaurus parameters object before running Thesaurus!");
		}
	}

	@Override
	protected Nothing process() {
		
		// break up decoys vs non-decoys by sequence
		HashMap<String, ArrayList<LibraryEntry>> entriesBySequence=new HashMap<>();
		for (LibraryEntry entry : super.entries) {
			String seq=entry.getPeptideSeq();
			ArrayList<LibraryEntry> list=entriesBySequence.get(seq);
			if (list==null) {
				list=new ArrayList<>();
				entriesBySequence.put(seq, list);
			}
			list.add(entry);
		}
		
		for (ArrayList<LibraryEntry> seedEntries : entriesBySequence.values()) {
			processPeptide(seedEntries);
		}
		return Nothing.NOTHING;
	}
	
	private void processPeptide(ArrayList<LibraryEntry> seedEntries) {
		LibraryEntry firstEntry=seedEntries.get(0);
		ArrayList<String> peptideModSeqs=PhosphoPermuter.getPermutations(firstEntry.getPeptideModSeq(), localizingModification, parameters.getAAConstants());

		byte precursorCharge=firstEntry.getPrecursorCharge();
		float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(firstEntry.getPeptideModSeq(), parameters.getAAConstants());

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
		for (String peptideModSeq : peptideModSeqs) {
			LocalizableForm form=getLocalizedForm(peptideModSeq, firstEntry.getPrecursorCharge(), entryMap, ionsByPeptide, seedEntries, parameters);
			if (form!=null) {
				allIsoforms.put(peptideModSeq, form);
				unlocalizedIsoforms.add(peptideModSeq);
			}
		}

		ArrayList<Spectrum>[] globalScans=getTargetSpectra(seedEntries);
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
		FragmentIonBlacklist takenIdentifiedIons=new FragmentIonBlacklist(parameters.getFragmentTolerance());
		
		AbstractScoringResult bestNonlocalizedResult=null;
		ModificationLocalizationData bestNonlocalizedData=null;
		boolean anyLocalized=false;
		
		HashMap<String, ScoredIndex> bestIndicies=new HashMap<>();
		TObjectIntHashMap<String> numberOfAttempts=new TObjectIntHashMap<>();
		
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
				
				LibraryEntry localizedEntry=localizedForm.localizedEntry;
				Ion[] allIons=localizedForm.allIons;
				

				for (int peak=0; peak<globalScans.length; peak++) {
					Float[] primary=globalPrimaryScores[peak].get(peptideModSeq);
					boolean[] alreadyConsidered=globalAlreadyConsideredScores[peak].get(peptideModSeq);

					if (primary==null) {
						primary=new Float[globalScans[peak].size()];
						alreadyConsidered=new boolean[primary.length];
						globalPrimaryScores[peak].put(peptideModSeq, primary);
						globalAlreadyConsideredScores[peak].put(peptideModSeq, alreadyConsidered);
					}
					ScoredIndex score=updateScores(globalScans[peak], localizedEntry, allIons, primary, alreadyConsidered, blacklistedScans, takenIdentifiedIons, parameters);
					
					ScoredIndex previousBest=bestIndicies.get(peptideModSeq);
					if (previousBest==null) {
						previousBest=score;
						bestIndicies.put(peptideModSeq, score);
					}

					if (score!=null&&unlocalizedIsoforms.contains(peptideModSeq)&&(bestIndex==null||score.x>bestIndex.x)&&(previousBest.x*MAXIMUM_DELTA_FROM_BEST_SCORE<score.x)) {
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
			numberOfAttempts.adjustOrPutValue(bestPeptideModSeq, 1, 1);
			ScoredIndex previousScoreIndex=bestIndicies.get(bestPeptideModSeq);
			if (previousScoreIndex!=null&&previousScoreIndex.x*MAXIMUM_DELTA_FROM_BEST_SCORE>=bestIndex.x) {
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
					
					FragmentIon[] localizingIons=getUniqueFragmentIons(entryMap.get(bestPeptideModSeq), entryMap.get(peptideModSeq), precursorCharge, parameters);
					
					Pair<FragmentScan, Float> localizedStripe = getBestLocalizationStripe(parameters, dutyCycle, localizer.getBackground(), bestForm.localizedEntry, localizingIons, stripeSubset);
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
				bestLocalizedStripe = getBestLocalizationStripe(parameters, dutyCycle, localizer.getBackground(), bestForm.localizedEntry, bestLocalizingIons, stripeSubset);
				if (bestLocalizedStripe==null||bestLocalizedStripe.x==null) continue;
			}
			
			AmbiguousPeptideModSeq ambiPeptideModSeq=AmbiguousPeptideModSeq.getUnambigous(bestPeptideModSeq, parameters.getLocalizingModification().get(), parameters.getAAConstants(), "");
			Triplet<ModificationLocalizationData, FragmentScan, Range> locData=generateLocalizationData(false, minimumScore, parameters, localizer.getModification(), bestForm.localizedEntry,
					ambiPeptideModSeq, bestLocalizingIons, bestForm.allIons, takenIdentifiedIons, stripeSubset, bestLocalizedStripe);
			if (locData==null) continue;
			
			// if localized, then keep and remove from localizedForms
			ModificationLocalizationData data=locData.x;
			FragmentScan apex=locData.y;
			Range peakRange=locData.z.addBuffer(dutyCycle);
			float score=((EncyclopediaScorer)scorer).score(bestForm.localizedEntry, apex, bestForm.allIons);
			boolean replaceBestNonLocalizedResult = bestNonlocalizedResult==null||(score>=bestNonlocalizedResult.getBestScore()&&bestNonlocalizedData.getLocalizationScore()<data.getLocalizationScore());
			//System.out.println(bestPeptideModSeq+" --> "+apex.getScanStartTime()/60f+" min, loc:"+locData.x.getLocalizationScore()+" vs "+bestLocalizedStripe.y+", score:"+score+" vs "+(bestNonlocalizedResult==null?"null":Float.toString(bestNonlocalizedResult.getBestScore()))+" ("+replaceBestNonLocalizedResult+")");
			if (replaceBestNonLocalizedResult||data.isLocalized()) {
				TFloatFloatHashMap scoreByRTMap=new TFloatFloatHashMap();
				Float[] primaryScoreArray=primaryScores.get(bestPeptideModSeq);
				for (int i=0; i<scans.size(); i++) {
					scoreByRTMap.put(scans.get(i).getScanStartTime(), primaryScoreArray[i]);
				}
				EValueCalculator calculator=new EValueCalculator(scoreByRTMap, 0f, 0.5f);
				float[] auxScoreArray=((EncyclopediaScorer)scorer).getAuxScorer().score(bestForm.localizedEntry, apex, predictedIsotopeDistribution, precursors);
	
				float evalue=calculator.getNegLnEValue(score);
				if (Float.isNaN(evalue)) {
					evalue=-1.0f;
				}

				AbstractScoringResult result=new PeptideScoringResult(bestForm.localizedEntry);

				float deltaPrecursorMass=scorer.getParentDeltaMassIndex()>=0?auxScoreArray[scorer.getParentDeltaMassIndex()]:0.0f;
				float deltaFragmentMass=scorer.getFragmentDeltaMassIndex()>=0?auxScoreArray[scorer.getFragmentDeltaMassIndex()]:0.0f;
				
				result.addStripe(score, General.concatenate(auxScoreArray, evalue, data.getLocalizationScore()), deltaPrecursorMass, deltaFragmentMass, apex);

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
				if (thesaurusParameters.isConsiderRearrangement()) {
					blacklistedScanRanges.remove(bestPeptideModSeq);
				}
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

			// TODO how do we know when to give up and just report a poor score?
			if (data.isLocalized()||numberOfAttempts.get(bestPeptideModSeq)>1) {
				// try up to 2x times
				unlocalizedIsoforms.remove(bestPeptideModSeq); // should we only do this if we can actually localize the peak?
			}
		}
		
		if (!anyLocalized&&bestNonlocalizedResult!=null) {
			resultsQueue.add(bestNonlocalizedResult);

			if (bestNonlocalizedResult.hasScoredResults()&&!bestNonlocalizedResult.getEntry().isDecoy()) {
				// don't bother logging decoys
				localizationQueue.add(bestNonlocalizedData);
			}
		}
	}

	public static ScoredIndex updateScores(ArrayList<Spectrum> scans, LibraryEntry localizedEntry, Ion[] allIons, Float[] primary, boolean[] alreadyConsidered, Collection<Range> blacklistedRanges, FragmentIonBlacklist takenIdentifiedIons, SearchParameters parameters) {
		float bestScore=-Float.MAX_VALUE;
		int bestIndex=0;
		for (int i=0; i<scans.size(); i++) {
			Spectrum stripe=scans.get(i);
			if (primary[i]==null) { 
				primary[i]=score(localizedEntry, stripe, allIons, takenIdentifiedIons, parameters);
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
	
	@SuppressWarnings("unchecked")
	private ArrayList<Spectrum>[] getTargetSpectra(ArrayList<LibraryEntry> seedEntries) {
		
		if (breadth==ScoringBreadthType.ENTIRE_RT_WINDOW) {
			return new ArrayList[] {PhosphoLocalizer.getScanSubsetFromStripes(-Float.MAX_VALUE, Float.MAX_VALUE, super.stripes)};

		} else if (breadth==ScoringBreadthType.UNCALIBRATED_20_PERCENT||breadth==ScoringBreadthType.RECALIBRATED_20_PERCENT) {
			float minTime=Float.MAX_VALUE;
			float maxTime=-Float.MAX_VALUE;
			for (LibraryEntry seed : seedEntries) {
				float rt;
				if (breadth==ScoringBreadthType.RECALIBRATED_20_PERCENT) {
					rt=seed.getScanStartTime();
				} else {
					rt=getBestScan(seed).getScanStartTime();
				}
				if (rt>maxTime) maxTime=rt;
				if (rt<minTime) minTime=rt;
			}
			float duration=(super.stripes.get(super.stripes.size()-1).getScanStartTime()-super.stripes.get(0).getScanStartTime())/20.0f; // 5% of entire gradient
			return new ArrayList[] {PhosphoLocalizer.getScanSubsetFromStripes(minTime-duration, maxTime+duration, super.stripes)};
			
		} else if (breadth==ScoringBreadthType.UNCALIBRATED_PEAK_WIDTH||breadth==ScoringBreadthType.RECALIBRATED_PEAK_WIDTH) {
			ArrayList<Spectrum>[] peaks=new ArrayList[seedEntries.size()];
			float duration=parameters.getExpectedPeakWidth()*0.75f; // 1.5x peak width total
			for (int i=0; i<peaks.length; i++) {
				float rt;
				if (breadth==ScoringBreadthType.RECALIBRATED_PEAK_WIDTH) {
					rt=getBestScan(seedEntries.get(i)).getScanStartTime();
				} else {
					rt=seedEntries.get(i).getScanStartTime();
				}
				peaks[i]=PhosphoLocalizer.getScanSubsetFromStripes(rt-duration, rt+duration, super.stripes);
			}
			return peaks;
		}

		throw new EncyclopediaException("Unexpected Thesaurus Scoring Breadth: "+breadth);
	}

	private Spectrum getBestScan(LibraryEntry seedEntry) {
		EncyclopediaScorer eScorer=(EncyclopediaScorer)scorer;
		FragmentationModel model=PeptideUtils.getPeptideModel(seedEntry.getPeptideModSeq(), parameters.getAAConstants());
		FragmentIon[] allIons=model.getPrimaryIonObjects(parameters.getFragType(), seedEntry.getPrecursorCharge(), true);
		float[] primary=new float[super.stripes.size()];
		for (int i=0; i<super.stripes.size(); i++) {
			Spectrum stripe=super.stripes.get(i);
			primary[i]=eScorer.score(seedEntry, stripe, allIons);
		}
		float bestScore=-Float.MAX_VALUE;
		Spectrum bestStripe=null;
		for (int i=0; i<primary.length; i++) {
			if (bestScore<primary[i]) {
				bestScore=primary[i];
				bestStripe=super.stripes.get(i);
			}
		}
		return bestStripe;
	}
	
	public static class LocalizableForm {
		public final LibraryEntry localizedEntry;
		public final FragmentIon[] allIons;
		public LocalizableForm(FragmentationModel localizedModel, LibraryEntry localizedEntry, SearchParameters parameters) {
			this.localizedEntry=localizedEntry;
			allIons=localizedModel.getPrimaryIonObjects(parameters.getFragType(), localizedEntry.getPrecursorCharge(), true);
		}
	}
	
	public static LocalizableForm getLocalizedForm(String peptideModSeq, byte charge, HashMap<String, FragmentationModel> modelMap, HashMap<String, double[]> ionsByPeptide, ArrayList<LibraryEntry> seedEntries, SearchParameters parameters) {
		FragmentationModel targetModel=modelMap.get(peptideModSeq);
		
		// first check if we have real data to cover this form
		for (LibraryEntry realEntry : seedEntries) {
			if (peptideModSeq.equals(realEntry.getAccuratePeptideModSeq(parameters.getAAConstants()))) {
				return new LocalizableForm(targetModel, realEntry, parameters);
			}
		}

		int bestNumMatching=-1;
		LibraryEntry bestRealEntry=null;
		
		double[] targetIons=ionsByPeptide.get(peptideModSeq);
		if (targetIons==null) {
			// This can happen when we're localizing around an n-term acetyl that's been rearranged incorrectly by reversing
			//Logger.errorLine("Missing target ions for: "+targetPeptideModSeq+", Found ions for: "+General.toString(ionsByPeptide.keySet())+", skipping form");
			return null;
		}
		for (LibraryEntry realEntry : seedEntries) {
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
			}
		}

		Pair<FragmentationModel, LibraryEntry> localizedForm=bestRealEntry.getEntryFromNewSequence(peptideModSeq, false, bestRealEntry.isDecoy(), bestRealEntry.isDecoy(), parameters);
		return new LocalizableForm(localizedForm.x, localizedForm.y, parameters);
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
	
	private static float score(LibraryEntry entry, Spectrum spectrum, Ion[] ions, FragmentIonBlacklist blacklistedIons, SearchParameters parameters) {
		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		MassTolerance libraryTolerance=parameters.getLibraryFragmentTolerance();

		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		int count=0;
		float dotProduct=0.0f;
		for (Ion targetIon : ions) {
			double target=targetIon.getMass();
			if (blacklistedIons.isBlacklisted(target, spectrum.getScanStartTime())) {
				continue;
			}
			
			int[] predictedIndicies=libraryTolerance.getIndicies(predictedMasses, target);
			float predictedIntensity=0.0f;
			for (int i=0; i<predictedIndicies.length; i++) {
				if (predictedIntensity<predictedIntensities[predictedIndicies[i]]) {
					predictedIntensity=predictedIntensities[predictedIndicies[i]];
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
				float peakScore=predictedIntensity*intensity;
				if (intensity>0.0f) {
					dotProduct+=peakScore;
					count++;
				}
			}
		}
		
		return Log.protectedLog10(dotProduct)+Log.logFactorial(count); // X!Tandem score
	}

	public static FragmentIon[] getUniqueFragmentIons(FragmentationModel target, FragmentationModel nextBest, byte precursorCharge, SearchParameters params) {
		HashSet<FragmentIon> ions=new HashSet<FragmentIon>(Arrays.asList(target.getPrimaryIonObjects(params.getFragType(), precursorCharge, false)));
		ions.removeAll(Arrays.asList(nextBest.getPrimaryIonObjects(params.getFragType(), precursorCharge, false)));

		FragmentIon[] ionArray=ions.toArray(new FragmentIon[ions.size()]);
		Arrays.sort(ionArray);
		return ionArray;
	}

	public static Triplet<ModificationLocalizationData, FragmentScan, Range> generateLocalizationData(
			boolean wasPreviouslyObservedSite, float minimumScore, SearchParameters parameters,
			PeptideModification mod, LibraryEntry localizedEntry, AmbiguousPeptideModSeq peptideModSeq,
			FragmentIon[] targetIons, FragmentIon[] allIons, FragmentIonBlacklist takenIdentifiedIons,
			ArrayList<Spectrum> stripeSubset, Pair<FragmentScan, Float> bestLocalizedStripe) {
		float bestLocalizationScore=bestLocalizedStripe.y;
		FragmentScan apex=bestLocalizedStripe.x;
		if (apex==null) {
			return null;
		}
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" --> ("+targetIons.length+") "+bestLocalizationScore+" localization score at "+(apex.getScanStartTime()/60f)+" min"); //FIXME 
		
		int targetNumFragments=Math.max(parameters.getMinNumOfQuantitativePeaks(), 3);

		boolean isLocalized=false;
		boolean isSiteSpecific=false;
		boolean isCompletelyAmbiguous=false;
		float localizationIntensity=0.0f;
		float totalIntensity=0.0f;
		float numIdentificationPeaks=0.0f;
		int numConsideredPeaks=0;
		int numberOfMods=PeptideUtils.getNumberOfMods(peptideModSeq.getPeptideModSeq(), mod.getNominalMass());
		ArrayList<FragmentIon> wellShapedIons=new ArrayList<FragmentIon>();

        Range peakRange=new Range(apex.getScanStartTime(), apex.getScanStartTime());
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score)\tNOT GOOD ENOUGH");
		if (bestLocalizationScore>=minimumScore||(bestLocalizationScore>0&&!wasPreviouslyObservedSite)) {
			// generate quant data from localizing ions only
			TransitionRefinementData quantData=PhosphoLocalizer.quantifyPeptide(parameters, peptideModSeq.getPeptideModSeq(), localizedEntry.getPrecursorCharge(), targetIons, apex.getScanStartTime(),
					stripeSubset, takenIdentifiedIons, Optional.ofNullable((float[]) null));
			//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score)\ta)"+(quantData!=null));
			if (quantData!=null) {
				peakRange=quantData.getRange();
				float[] intensities=quantData.getIntegrationArray();
				float[] correlations=quantData.getCorrelationArray();
				FragmentIon[] consideredIons=(FragmentIon[])quantData.getFragmentMassArray();
				for (int i=0; i<consideredIons.length; i++) {
					wellShapedIons.add(consideredIons[i]);
					localizationIntensity+=intensities[i];
				}

				// calculate quant data for all ions
				float[] medianChromatogram=quantData.getMedianChromatogram();
				TransitionRefinementData allQuantData=PhosphoLocalizer.quantifyPeptide(parameters, peptideModSeq.getPeptideModSeq(), localizedEntry.getPrecursorCharge(), allIons, apex.getScanStartTime(),
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
							matchingPeaks.add(allIons[i].getMass());
						}
						//System.out.println(localizedEntry.getPeptideModSeq()+", "+numIdentificationPeaks+") "+targetIonSet.contains(allIons[i])+", "+allIons[i].toString()+", "+allIons[i].mass+"\t"+intensities[i]+"\t"+correlations[i]);
					}
					// recalculate localization scores using only the peaks are remotely ok
					// FIXME
//					double[] targetIonsMasses=FragmentIon.getMasses(targetIons);
//					float[] frequencies=localizer.getBackground().getFrequencies(targetIonsMasses, localizedEntry.getPrecursorMZ(), parameters.getFragmentTolerance());
//					bestLocalizationScore=PhosphoLocalizer.score(parameters, targetIonsMasses, targetIons, frequencies, matchingPeaks.toArray());

					//System.out.println(localizedEntry.getPeptideModSeq()+" Final: "+bestLocalizationScore);
					
					isCompletelyAmbiguous=AmbiguousPeptideModSeq.isCompletelyAmbiguous(peptideModSeq, mod);
					isLocalized=bestLocalizationScore>=minimumScore&&wellShapedIons.size()>0&&numIdentificationPeaks>targetNumFragments&&!isCompletelyAmbiguous;
					isSiteSpecific=isLocalized&&AmbiguousPeptideModSeq.isSiteSpecific(peptideModSeq, mod);
					//System.out.println("\tLocalized "+isLocalized+"/"+isSiteSpecific+" for "+peptideModSeq.getPeptideAnnotation()+" ("+bestLocalizationScore+" score vs "+minimumScore+"minimum, "+numIdentificationPeaks+"/"+correlations.length+" peaks)"); // FIXME
				}
			}
		}
		//System.out.println("CASiL: "+peptideModSeq.getPeptideModSeq()+" ("+targetIons.length+"/"+allIons.length+"/"+numConsideredPeaks+")\t"+bestLocalizationScore+","+isCompletelyAmbiguous+","+isLocalized+","+isSiteSpecific+","+numIdentificationPeaks);
		
		ModificationLocalizationData modData=new ModificationLocalizationData(peptideModSeq, apex.getScanStartTime(), bestLocalizationScore, numIdentificationPeaks, numberOfMods, isSiteSpecific, isLocalized, isCompletelyAmbiguous, wellShapedIons.toArray(new FragmentIon[wellShapedIons.size()]), localizationIntensity, totalIntensity);
		return new Triplet<ModificationLocalizationData, FragmentScan, Range>(modData, apex, peakRange);
	}

	public static Pair<FragmentScan, Float> getBestLocalizationStripe(SearchParameters parameters, float dutyCycle,
			BackgroundFrequencyInterface background, LibraryEntry localizedEntry, FragmentIon[] targetIons,
			ArrayList<Spectrum> stripeSubset) {
		double[] targetIonsMasses=FragmentIon.getMasses(targetIons);
		float[] frequencies=background.getFrequencies(targetIonsMasses, localizedEntry.getPrecursorMZ(), parameters.getFragmentTolerance());
		
		//System.out.println("\t"+peptideModSeq.getPeptideAnnotation()+" --> "+General.toString(targetIonsMasses));

		float[] negLogProbsSiteSpecific=new float[stripeSubset.size()];
		for (int k=0; k<stripeSubset.size(); k++) {
			negLogProbsSiteSpecific[k]=PhosphoLocalizer.score(parameters, targetIonsMasses, targetIons, frequencies, stripeSubset.get(k), true);
		}
		float bestScore=-Float.MAX_VALUE;
		FragmentScan apexStripe=null;
		negLogProbsSiteSpecific=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, Math.round(parameters.getExpectedPeakWidth()/dutyCycle));//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);
		for (int k=0; k<stripeSubset.size(); k++) {
			if (bestScore<negLogProbsSiteSpecific[k]) {
				bestScore=negLogProbsSiteSpecific[k];
				apexStripe=(FragmentScan)stripeSubset.get(k);
			}
		}
		
		Pair<FragmentScan, Float> bestLocalizedStripe=new Pair<FragmentScan, Float>(apexStripe, bestScore);
		return bestLocalizedStripe;
	}
}
