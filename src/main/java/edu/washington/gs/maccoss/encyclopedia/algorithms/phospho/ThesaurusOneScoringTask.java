package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.IsotopicDistributionCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class ThesaurusOneScoringTask extends AbstractLibraryScoringTask {
	
	private final PhosphoLocalizer localizer;
	private final float dutyCycle;
	private final ScoringBreadthType breadth;
	private final PeptideModification localizingModification;
	private final BlockingQueue<ModificationLocalizationData> localizationQueue;
	private final float minimumScore;
	
	public ThesaurusOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, 
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
			FragmentationModel model=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
			entryMap.put(peptideModSeq, model);
			ionsByPeptide.put(peptideModSeq, model.getPrimaryIons(parameters.getFragType(), firstEntry.getPrecursorCharge()));
		}

		HashSet<String> unlocalizedIsoforms=new HashSet<>();
		HashMap<String, LocalizableForm> allIsoforms=new HashMap<>();
		for (String peptideModSeq : peptideModSeqs) {
			LocalizableForm form=getLocalizedForm(peptideModSeq, firstEntry.getPrecursorCharge(), entryMap, ionsByPeptide, seedEntries, parameters);
			if (form!=null) {
				allIsoforms.put(peptideModSeq, form);
				unlocalizedIsoforms.add(peptideModSeq);
			}
		}

		HashMap<String, Float[]> primaryScores=new HashMap<>(); // scores are re-used unless they fall in the RT range of a previously localized form
		ArrayList<Spectrum> scans=getTargetSpectra(seedEntries);
		FragmentIonBlacklist takenIdentifiedIons=new FragmentIonBlacklist(parameters.getFragmentTolerance());
		
		PeptideScoringResult bestNonlocalizedResult=null;
		ModificationLocalizationData bestNonlocalizedData=null;
		boolean anyLocalized=false;
		
		while (unlocalizedIsoforms.size()>0) {
			String bestPeptideModSeq=null;
			ScoredIndex bestIndex=null;
			LocalizableForm bestForm=null;
			
			// get the highest scoring RT point for all peptide sequences
			for (Entry<String, LocalizableForm> entry : allIsoforms.entrySet()) {
				String peptideModSeq=entry.getKey();
				LocalizableForm localizedForm=entry.getValue();
				
				LibraryEntry localizedEntry=localizedForm.localizedEntry;
				FragmentIon[] allIons=localizedForm.allIons;
				
				Float[] primary=primaryScores.get(peptideModSeq);
				if (primary==null) {
					primary=new Float[scans.size()];
					primaryScores.put(peptideModSeq, primary);
				}
				ScoredIndex score=updateScores(scans, localizedEntry, allIons, primary, takenIdentifiedIons);

				if (unlocalizedIsoforms.contains(peptideModSeq)&&(bestIndex==null||score.x>bestIndex.x)) {
					bestIndex=score;
					bestPeptideModSeq=peptideModSeq;
					bestForm=localizedForm;
				}
			}
			System.out.println("CHECK: "+bestPeptideModSeq+"\t"+bestIndex.x+"\t"+(scans.get(bestIndex.y).getScanStartTime()/60f)+" (total scans: "+scans.size()+")"); //FIXME
			ArrayList<XYTraceInterface> traces=new ArrayList<>();
			float[] rts=new float[scans.size()];
			for (int i = 0; i < rts.length; i++) {
				rts[i]=scans.get(i).getScanStartTime()/60f;
			}
			for (Entry<String, Float[]> scores : primaryScores.entrySet()) {
				traces.add(new XYTrace(rts, General.toFloatArray(scores.getValue()), GraphType.boldline, scores.getKey()));
			}
			Charter.launchChart("RT", "Score", true, traces.toArray(new XYTraceInterface[traces.size()]));
			
			// get next best match at that RT
			float nextBestScore=-Float.MAX_VALUE;
			String nextBestPeptideModSeq=null;
			for (Entry<String, Float[]> entry : primaryScores.entrySet()) {
				String peptideModSeq=entry.getKey();
				if (peptideModSeq!=bestPeptideModSeq) {
					float score=entry.getValue()[bestIndex.y];
					if (score!=THIS_PEPTIDE_IS_NOT_HERE&&score>nextBestScore) {
						nextBestPeptideModSeq=peptideModSeq;
						nextBestScore=score;
					}
				}
			}
			
			// check localization ions versus that sequence
			float apexRT=scans.get(bestIndex.y).getScanStartTime();
			// use stripes here in case we're on the border
			ArrayList<Spectrum> stripeSubset=PhosphoLocalizer.getScanSubsetFromStripes(apexRT-parameters.getExpectedPeakWidth(), apexRT+parameters.getExpectedPeakWidth(), stripes); 
			
			FragmentIon[] localizingIons;
			if (nextBestPeptideModSeq!=null) {
				localizingIons=getUniqueFragmentIons(entryMap.get(bestPeptideModSeq), entryMap.get(nextBestPeptideModSeq), precursorCharge, parameters);
			} else {
				localizingIons=bestForm.allIons;
			}
			AmbiguousPeptideModSeq ambiPeptideModSeq=AmbiguousPeptideModSeq.getUnambigous(bestPeptideModSeq, parameters.getLocalizingModification().get(), parameters.getAAConstants());
			Triplet<ModificationLocalizationData, Stripe, Range> locData=CASiLOneScoringTask.calculateLocalizationScoring(false, minimumScore, parameters, dutyCycle, localizer, bestForm.localizedEntry, ambiPeptideModSeq, localizingIons, bestForm.allIons, takenIdentifiedIons, stripeSubset);
			
			// if localized, then keep and remove from localizedForms
			ModificationLocalizationData data=locData.x;
			Stripe apex=locData.y;
			Range peakRange=locData.z;
			
			boolean wasLocalized=bestNonlocalizedResult==null||data.isLocalized();
			if (wasLocalized) {
				
				TFloatFloatHashMap scoreByRTMap=new TFloatFloatHashMap();
				Float[] primaryScoreArray=primaryScores.get(bestPeptideModSeq);
				for (int i=0; i<scans.size(); i++) {
					scoreByRTMap.put(scans.get(i).getScanStartTime(), primaryScoreArray[i]);
				}
				EValueCalculator calculator=new EValueCalculator(scoreByRTMap);
				float[] auxScoreArray=((EncyclopediaScorer)scorer).getAuxScorer().score(bestForm.localizedEntry, apex, predictedIsotopeDistribution, precursors);
	
				float score=((EncyclopediaScorer)scorer).score(bestForm.localizedEntry, apex, bestForm.allIons);
				float evalue=calculator.getNegLog10EValue(score);
				if (Float.isNaN(evalue)) {
					evalue=-1.0f;
				}

				PeptideScoringResult result=new PeptideScoringResult(bestForm.localizedEntry);
				result.addStripe(score, General.concatenate(auxScoreArray, evalue, data.getLocalizationScore()), apex);

				if (bestNonlocalizedResult==null) {
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

			// FIXME how do we know when to give up and just report a poor score? 
			unlocalizedIsoforms.remove(bestPeptideModSeq); // should we only do this if we can actually localize the peak?
			
			for (FragmentIon target : localizingIons) {
				takenIdentifiedIons.addIonToBlacklist(target.mass, peakRange);
			}
			// null out scores from taken ions
			for (int i=0; i<scans.size(); i++) {
				if (peakRange.contains(scans.get(i).getScanStartTime())) {
					for (Entry<String, Float[]> entry : primaryScores.entrySet()) {
						if ((!wasLocalized)||entry.getKey()!=bestPeptideModSeq) {
							entry.getValue()[i]=null;
						}
					}
				} else {
					if (wasLocalized) {
						// we found it elsewhere so it wasn't here
						primaryScores.get(bestPeptideModSeq)[i]=THIS_PEPTIDE_IS_NOT_HERE;
					}
				}
			}
		}
		
		if (!anyLocalized&&bestNonlocalizedResult!=null) {
			resultsQueue.add(bestNonlocalizedResult);

			if (!bestNonlocalizedResult.getEntry().isDecoy()) {
				// don't bother logging decoys
				localizationQueue.add(bestNonlocalizedData);
			}
		}
	}
	private final Float THIS_PEPTIDE_IS_NOT_HERE=-1f;

	private ScoredIndex updateScores(ArrayList<Spectrum> scans, LibraryEntry localizedEntry, FragmentIon[] allIons, Float[] primary, FragmentIonBlacklist takenIdentifiedIons) {
		ScoredIndex bestIndex=null;
		for (int i=0; i<scans.size(); i++) {
			Spectrum stripe=scans.get(i);
			if (primary[i]==null) { 
				primary[i]=score(localizedEntry, stripe, allIons, takenIdentifiedIons, parameters);
			}
			if (bestIndex==null||bestIndex.x<primary[i]) {
				bestIndex=new ScoredIndex(primary[i], i);
			}
		}
		return bestIndex;
	}
	
	private ArrayList<Spectrum> getTargetSpectra(ArrayList<LibraryEntry> seedEntries) {
		
		if (breadth==ScoringBreadthType.ENTIRE_RT_WINDOW) {
			return PhosphoLocalizer.getScanSubsetFromStripes(-Float.MAX_VALUE, Float.MAX_VALUE, super.stripes);

		} else if (breadth==ScoringBreadthType.UNCALIBRATED_20_PERCENT||breadth==ScoringBreadthType.UNCALIBRATED_PEAK_WIDTH) {
			float minTime=Float.MAX_VALUE;
			float maxTime=-Float.MAX_VALUE;
			for (LibraryEntry seed : seedEntries) {
				if (seed.getScanStartTime()>maxTime) maxTime=seed.getScanStartTime();
				if (seed.getScanStartTime()<minTime) minTime=seed.getScanStartTime();
			}
			
			float duration;
			if (breadth==ScoringBreadthType.UNCALIBRATED_20_PERCENT) {
				duration=(super.stripes.get(super.stripes.size()-1).getScanStartTime()-super.stripes.get(0).getScanStartTime())/20.0f;
			} else { //if (breadth==ScoringBreadthType.UNCALIBRATED_PEAK_WIDTH) 
				duration=parameters.getExpectedPeakWidth();
			}
			return PhosphoLocalizer.getScanSubsetFromStripes(minTime-duration, maxTime+duration, super.stripes);
			
		} else if (breadth==ScoringBreadthType.RECALIBRATED_20_PERCENT||breadth==ScoringBreadthType.RECALIBRATED_PEAK_WIDTH) {
			EncyclopediaScorer eScorer=(EncyclopediaScorer)scorer;
			
			float minTime=Float.MAX_VALUE;
			float maxTime=-Float.MAX_VALUE;
			for (LibraryEntry seedEntry : seedEntries) {
				FragmentationModel model=new FragmentationModel(seedEntry.getPeptideModSeq(), parameters.getAAConstants());
				FragmentIon[] allIons=model.getPrimaryIonObjects(parameters.getFragType(), seedEntry.getPrecursorCharge());
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
				if (bestStripe.getScanStartTime()>maxTime) maxTime=bestStripe.getScanStartTime();
				if (bestStripe.getScanStartTime()<minTime) minTime=bestStripe.getScanStartTime();
			}
			
			float duration;
			if (breadth==ScoringBreadthType.RECALIBRATED_20_PERCENT) {
				duration=(super.stripes.get(super.stripes.size()-1).getScanStartTime()-super.stripes.get(0).getScanStartTime())/20.0f;
			} else { // if (breadth==ScoringBreadthType.RECALIBRATED_PEAK_WIDTH)
				duration=parameters.getExpectedPeakWidth();
			}
			
			return PhosphoLocalizer.getScanSubsetFromStripes(minTime-duration, maxTime+duration, super.stripes);
		}

		throw new EncyclopediaException("Unexpected CASiL Scoring Breadth: "+breadth);
	}
	
	static class LocalizableForm {
		private final LibraryEntry localizedEntry;
		private final FragmentIon[] allIons;
		public LocalizableForm(FragmentationModel localizedModel, LibraryEntry localizedEntry, SearchParameters parameters) {
			this.localizedEntry=localizedEntry;
			allIons=localizedModel.getPrimaryIonObjects(parameters.getFragType(), localizedEntry.getPrecursorCharge());
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

		Pair<FragmentationModel, LibraryEntry> localizedForm=bestRealEntry.getEntryFromNewSequence(peptideModSeq, bestRealEntry.getAccessions(), bestRealEntry.isDecoy(), parameters);
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
	
	private static float score(LibraryEntry entry, Spectrum spectrum, FragmentIon[] ions, FragmentIonBlacklist blacklistedIons, SearchParameters parameters) {
		MassTolerance acquiredTolerance=parameters.getFragmentTolerance();
		MassTolerance libraryTolerance=parameters.getLibraryFragmentTolerance();

		double[] predictedMasses=entry.getMassArray();
		float[] predictedIntensities=entry.getIntensityArray();
		
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();

		int count=0;
		float dotProduct=0.0f;
		for (FragmentIon targetIon : ions) {
			double target=targetIon.mass;
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
}
