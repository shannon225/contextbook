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
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideScoringResult;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Nothing;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredIndex;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.set.hash.TIntHashSet;

public class CAPSiLOneScoringTask extends AbstractLibraryScoringTask {
	private final PhosphoLocalizer localizer;
	private final float dutyCycle;
	private final int targetNumFragments;
	
	public CAPSiLOneScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<Stripe> stripes, float dutyCycle, PrecursorScanMap precursors, 
			PhosphoLocalizer localizer, BlockingQueue<PeptideScoringResult> resultsQueue, SearchParameters parameters) {
		super(scorer, entries, stripes, precursors, resultsQueue, parameters);
		this.dutyCycle=dutyCycle;
		this.localizer=localizer;
		targetNumFragments=Math.max(parameters.getMinNumOfQuantitativePeaks(), 3);
		
	}
	private static final int peaksKept=-1; // only keep the first peak

	@Override
	protected Nothing process() {
		EncyclopediaScorer eScorer=(EncyclopediaScorer)scorer;
		int movingAverageLength=Math.round(parameters.getExpectedPeakWidth()/dutyCycle);
		for (LibraryEntry seedEntry : super.entries) {
			float[] predictedIsotopeDistribution=IsotopicDistributionCalculator.getIsotopeDistribution(seedEntry.getPeptideModSeq(), parameters.getAAConstants());
			
			ArrayList<String> peptideModSeqs=PhosphoPermuter.getPermutations(seedEntry.getPeptideModSeq(), parameters.getAAConstants());

			HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
			for (String peptideModSeq : peptideModSeqs) {
				FragmentationModel model=new FragmentationModel(peptideModSeq, parameters.getAAConstants());
				entryMap.put(peptideModSeq, model);
			}

			TIntHashSet takenScans=new TIntHashSet();
			while (peptideModSeqs.size()!=0) {
				boolean breakBatch=false;

				ArrayList<Pair<AmbiguousPeptideModSeq, FragmentIon[]>> batch=new ArrayList<Pair<AmbiguousPeptideModSeq,FragmentIon[]>>();
				AmbiguousPeptideModSeq leftAmbiguity=AmbiguousPeptideModSeq.getLeftAmbiguity(peptideModSeqs.remove(0), AmbiguousPeptideModSeq.modifiableAAs, parameters.getAAConstants());
				FragmentIon[] leftTargets=PhosphoLocalizer.getUniqueFragmentIons(leftAmbiguity.getPeptideModSeq(), seedEntry.getPrecursorCharge(), entryMap, parameters);
				batch.add(new Pair<AmbiguousPeptideModSeq, FragmentIon[]>(leftAmbiguity, leftTargets));
				if (peptideModSeqs.size()>0) {
					AmbiguousPeptideModSeq rightAmbiguity=AmbiguousPeptideModSeq.getRightAmbiguity(peptideModSeqs.remove(peptideModSeqs.size()-1), AmbiguousPeptideModSeq.modifiableAAs, parameters.getAAConstants());
					FragmentIon[] rightTargets=PhosphoLocalizer.getUniqueFragmentIons(rightAmbiguity.getPeptideModSeq(), seedEntry.getPrecursorCharge(), entryMap, parameters);
					batch.add(new Pair<AmbiguousPeptideModSeq, FragmentIon[]>(rightAmbiguity, rightTargets));
				}
								
				for (Pair<AmbiguousPeptideModSeq, FragmentIon[]> pair : batch) {
					AmbiguousPeptideModSeq peptideModSeq=pair.x;
					FragmentIon[] targetIons=pair.y;
					
					entryMap.remove(peptideModSeq); // allows searching beyond this mod (only if we localize it to a RT)
									
					Pair<FragmentationModel, LibraryEntry> localizedForm=seedEntry.getEntryFromNewSequence(peptideModSeq.getPeptideModSeq(), seedEntry.getAccessions(), seedEntry.isDecoy(), parameters);
					FragmentationModel localizedModel=localizedForm.x;
					LibraryEntry localizedEntry=localizedForm.y;

					AuxillaryPSMScorer auxScorer=eScorer.getAuxScorer().getEntryOptimizedScorer(localizedEntry);
					FragmentIon[] allIons=localizedModel.getPrimaryIonObjects(parameters.getFragType(), localizedEntry.getPrecursorCharge());

					float[] primary=new float[super.stripes.size()];
					for (int i=0; i<super.stripes.size(); i++) {
						Stripe stripe=super.stripes.get(i);
						primary[i]=eScorer.score(localizedEntry, stripe, allIons);
					}
					
					float[] averagePrimary=gaussianCenteredAverage(primary, movingAverageLength);

					TFloatFloatHashMap scoreByRTMap=new TFloatFloatHashMap();
					ArrayList<ScoredIndex> goodStripes=new ArrayList<ScoredIndex>();
					for (int i=0; i<averagePrimary.length; i++) {
						goodStripes.add(new ScoredIndex(primary[i], i));
						scoreByRTMap.put(super.stripes.get(i).getScanStartTime(), primary[i]);
					}
					Collections.sort(goodStripes);

					EValueCalculator calculator=new EValueCalculator(scoreByRTMap);

					//System.out.println("Considering "+localizedEntry.getPeptideModSeq()+"\t from seed: "+seedEntry.getPeptideModSeq()); // FIXME
					PeptideScoringResult result=new PeptideScoringResult(localizedEntry);
					// find the best stripe and localize near it
					int identifiedPeaks=0;
					for (int i=goodStripes.size()-1; i>=0; i--) {
						int index=goodStripes.get(i).y;
						if (takenScans.contains(index)) {
							continue;
							
						} else {
							Stripe stripe=super.stripes.get(index);
							
							ArrayList<Spectrum> stripeSubset=PhosphoLocalizer.getScanSubsetFromStripes(stripe.getScanStartTime()-parameters.getExpectedPeakWidth(), stripe.getScanStartTime()+parameters.getExpectedPeakWidth(), super.stripes);
														
							double[] targetIonsMasses=FragmentIon.getMasses(targetIons);
							float[] frequencies=localizer.getBackground().getFrequencies(targetIonsMasses, seedEntry.getPrecursorMZ(), parameters.getFragmentTolerance());
							
							float bestLocalizationScore=-Float.MAX_VALUE;
							Stripe apex=stripe;
							float[] negLogProbsSiteSpecific=new float[stripeSubset.size()];
							for (int k=0; k<stripeSubset.size(); k++) {
								negLogProbsSiteSpecific[k]=PhosphoLocalizer.score(parameters, targetIonsMasses, targetIons, frequencies, stripeSubset.get(k), true);
								if (bestLocalizationScore<negLogProbsSiteSpecific[k]) {
									bestLocalizationScore=negLogProbsSiteSpecific[k];
									apex=(Stripe)stripeSubset.get(k);
								}
							}
							//System.out.println("\t("+i+")"+peptideModSeq.getPeptideAnnotation()+" --> "+bestScore+" localization score at "+(apex.getScanStartTime()/60f)+" min"); //FIXME 

							boolean wasLocalized=false;
							if (bestLocalizationScore>=PhosphoLocalizer.MINIMUM_SCORE) {
								// generate quant data from localizing ions only
								TransitionRefinementData quantData=localizer.quantifyPeptide(peptideModSeq.getPeptideModSeq(), seedEntry.getPrecursorCharge(), targetIons, apex.getScanStartTime(),
										stripeSubset, Optional.ofNullable((float[]) null));
								if (quantData!=null) {

									// calculate quant data for all ions
									float[] medianChromatogram=quantData.getMedianChromatogram();
									TransitionRefinementData allQuantData=localizer.quantifyPeptide(peptideModSeq.getPeptideModSeq(), seedEntry.getPrecursorCharge(), allIons, apex.getScanStartTime(),
											stripeSubset, Optional.of(medianChromatogram));
									if (allQuantData!=null) {

										int numIdentificationPeaks=0;
										float[] correlations=allQuantData.getCorrelationArray();
										for (int k=0; k<correlations.length; k++) {
											if (correlations[k]>=TransitionRefiner.identificationCorrelationThreshold) {
												numIdentificationPeaks++;
											}
										}
										wasLocalized=numIdentificationPeaks>=targetNumFragments&&AmbiguousPeptideModSeq.isLocalized(peptideModSeq, AmbiguousPeptideModSeq.modifiableAAs);
										//System.out.println("\tLocalized "+wasLocalized+" for "+peptideModSeq.getPeptideAnnotation()+" ("+bestScore+" score, "+numIdentificationPeaks+"/"+correlations.length+" peaks)"); // FIXME

										if (wasLocalized) {
											float[] auxScoreArray=auxScorer.score(localizedEntry, apex, predictedIsotopeDistribution, precursors);

											float score=eScorer.score(localizedEntry, apex, allIons);
											float evalue=calculator.getNegLog10EValue(score);
											if (Float.isNaN(evalue)) {
												evalue=-1.0f;
											}

											result.addStripe(score, General.concatenate(auxScoreArray, evalue, bestLocalizationScore), apex);
											resultsQueue.add(result);
										}
									}
								}
							}
							if (!wasLocalized) {
								breakBatch=true;
							}
							
							// block out a half a peakWidth window
							int removedWindow=Math.max(1, movingAverageLength/4);
							int lowerWindow=index-removedWindow;
							int upperWindow=index+removedWindow;
							for (int j=lowerWindow; j<=upperWindow; j++) {
								takenScans.add(j);
							}
							
							//System.out.println("\t\tcheck "+identifiedPeaks+" > "+peaksKept); // FIXME
							if (identifiedPeaks>peaksKept) {
								// keep N+1 peaks
								break;
							}
							identifiedPeaks++;
						}
					}
				}
				if (breakBatch) {
					// wasn't able to localize a site, so can't keep going
					break;
				}
			}
		}
		return Nothing.NOTHING;
	}
	
}
