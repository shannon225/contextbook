package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ModificationLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PhosphoLocalizer {
	private final float minimumScore;
	private final StripeFileInterface diaFile;
	private final SearchParameters params;
	private final BackgroundFrequencyInterface background;
	private final float gradientLength;
	private final PeptideModification modification;

	public PhosphoLocalizer(StripeFileInterface diaFile, PeptideModification localizingModification, SearchParameters params) throws IOException,DataFormatException,SQLException {
		this.diaFile=diaFile;
		this.params=params;
		this.modification=localizingModification;
		this.minimumScore=-Log.log10(params.getPercolatorThreshold());
		this.background=BackgroundFrequencyCalculator.generateBackground(diaFile);
		this.gradientLength=diaFile.getGradientLength();
	}

	public PhosphoLocalizer(StripeFileInterface diaFile, PeptideModification localizingModification, BackgroundFrequencyInterface background, SearchParameters params) throws IOException,DataFormatException,SQLException {
		this.diaFile=diaFile;
		this.modification=localizingModification;
		this.params=params;
		this.minimumScore=-Log.log10(params.getPercolatorThreshold());
		this.background=background;
		this.gradientLength=diaFile.getGradientLength();
	}

	public PhosphoLocalizer(StripeFileInterface diaFile, PeptideModification localizingModification, LibraryInterface searchedLibrary, SearchParameters params) throws IOException,DataFormatException,SQLException {
		this.diaFile=diaFile;
		this.modification=localizingModification;
		this.params=params;
		this.minimumScore=-Log.log10(params.getPercolatorThreshold());
		this.background=BackgroundFrequencyCalculator.generateBackground(diaFile);
		this.gradientLength=diaFile.getGradientLength();
	}
	
	public PeptideModification getModification() {
		return modification;
	}
	
	public BackgroundFrequencyInterface getBackground() {
		return background;
	}

	public Optional<PhosphoLocalizationData> runDIAPhosphoLocalization(PSMData psmdata, ArrayList<Stripe> stripes, boolean tryAllPermutations, boolean buildOutGUIData) {
		ArrayList<Spectrum> spectra=new ArrayList<Spectrum>();
		for (Stripe stripe : stripes) {
			spectra.add(stripe);
		}
		return runPhosphoLocalization(psmdata, spectra, tryAllPermutations, buildOutGUIData);
	}

	public Optional<PhosphoLocalizationData> runPhosphoLocalization(PSMData psmdata, ArrayList<Spectrum> stripes, boolean tryAllPermutations, boolean buildOutGUIData) {
		ArrayList<String> permutations;
		if (tryAllPermutations) {
			permutations=PhosphoPermuter.getPermutations(psmdata.getPeptideModSeq(), modification, params.getAAConstants());
		} else {
			permutations=new ArrayList<String>();
			permutations.add(psmdata.getPeptideModSeq());
		}
		if (permutations.size()==1) {
			//System.out.println("single\t"+psmdata.getPeptideModSeq()+"\t1\t1\t0\t1000");
			return Optional.empty();
		} else {
			PhosphoLocalizationData multiple=extractPhosphoForms(psmdata.getPeptideModSeq(), psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), permutations, psmdata.getRetentionTime(), stripes, buildOutGUIData);
			return Optional.of(multiple);
		}
	}

	/**
	 * 
	 * @param precursorMZ
	 * @param precursorCharge
	 * @param peptideModSeqs sequences are in phospho order left to right
	 * @param retentionTime
	 * @param allScansInStripe
	 * @return
	 */
	PhosphoLocalizationData extractPhosphoForms(String originalPeptideModSeq, double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, float retentionTime, ArrayList<Spectrum> allScansInStripe, boolean buildOutGUIData) {
		float duration=gradientLength/20.0f;
		ArrayList<Spectrum> stripes=getScanSubset(retentionTime-duration, retentionTime+duration, allScansInStripe);
		
		return extractPhosphoFormsFromLimitedScans(precursorMZ, precursorCharge, peptideModSeqs, stripes, buildOutGUIData);
	}

	PhosphoLocalizationData extractPhosphoFormsFromStripes(String originalPeptideModSeq, double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, float retentionTime, ArrayList<Stripe> allScansInStripe, boolean buildOutGUIData) {
		float duration=gradientLength/20.0f;
		ArrayList<Spectrum> stripes=getScanSubsetFromStripes(retentionTime-duration, retentionTime+duration, allScansInStripe);
		
		return extractPhosphoFormsFromLimitedScans(precursorMZ, precursorCharge, peptideModSeqs, stripes, buildOutGUIData);
	}

	public PhosphoLocalizationData extractPhosphoFormsFromLimitedScans(double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, ArrayList<Spectrum> stripes, boolean buildOutGUIData) {
		float dutyCycle=1.0f;
		for (Entry<Range, Float> entry : diaFile.getRanges().entrySet()) {
			if (entry.getKey().contains((float)precursorMZ)) {
				dutyCycle=entry.getValue();
				break;
			}
		}
		
		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
		for (String peptideModSeq : peptideModSeqs) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			entryMap.put(peptideModSeq, model);
		}
		
		// Iteratively subtract out peaks in either direction.
		// For SGS[+80]VS[+80]NYR, process:
		// (S[+80])GSVSNYR, (SGS[+80])VSNYR, (SGSVS[+80])NYR (but drop the last)
		// then from right to left:
		// (S[+80]GSVS)NYR, SG(S[+80]VS)NYR, SGSV(S[+80])NYR (but drop the first)
		//
		// For multiply phosphorylated SGS[+80]VS[+80]NYR, process:
		// (S[+80])G(S[+80])VSNYR, (S[+80])G(SVS[+80])NYR, (SGS[+80])V(S[+80])NYR (but drop the last)
		// then from right to left:
		// (S[+80])G(S[+80]VS)NYR, (S[+80]GS)V(S[+80])NYR, SG(S[+80])V(S[+80])NYR (but drop the first)
		
		ArrayList<Pair<AmbiguousPeptideModSeq, FragmentIon[]>> targetPeptidesLeft=new ArrayList<Pair<AmbiguousPeptideModSeq,FragmentIon[]>>(); 
		// go left to right, drop the last
		for (int i=0; i<peptideModSeqs.size()-1; i++) {
			String targetPeptide=peptideModSeqs.get(i);
			AmbiguousPeptideModSeq targetPeptideName=AmbiguousPeptideModSeq.getLeftAmbiguity(targetPeptide, modification, params.getAAConstants());

			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of
			// remaining sequences to the right
			for (int j=peptideModSeqs.size()-1; j>=i; j--) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			FragmentIon[] targets=PhosphoLocalizer.getUniqueFragmentIons(targetPeptide, precursorCharge, modelBatch, params);
			targetPeptidesLeft.add(new Pair<AmbiguousPeptideModSeq, FragmentIon[]>(targetPeptideName, targets));
		}
		ArrayList<Pair<AmbiguousPeptideModSeq, FragmentIon[]>> targetPeptidesRight=new ArrayList<Pair<AmbiguousPeptideModSeq,FragmentIon[]>>(); 
		// go right to left, drop the first
		for (int i=peptideModSeqs.size()-1; i>=1; i--) {
			String targetPeptide=peptideModSeqs.get(i);
			AmbiguousPeptideModSeq targetPeptideName=AmbiguousPeptideModSeq.getRightAmbiguity(targetPeptide, modification, params.getAAConstants());
			
			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of remaining sequences to the left
			for (int j=0; j<=i; j++) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			FragmentIon[] targets=PhosphoLocalizer.getUniqueFragmentIons(targetPeptide, precursorCharge, modelBatch, params);
			targetPeptidesRight.add(new Pair<AmbiguousPeptideModSeq, FragmentIon[]>(targetPeptideName, targets));
		}
		
		// interlay the peptides so we look at the most localizing first
		ArrayList<Pair<AmbiguousPeptideModSeq, FragmentIon[]>> targetPeptides=new ArrayList<Pair<AmbiguousPeptideModSeq,FragmentIon[]>>();
		for (int i=0; i<targetPeptidesLeft.size(); i++) {
			targetPeptides.add(targetPeptidesLeft.get(i));
			targetPeptides.add(targetPeptidesRight.get(i));
		}
		
		for (Pair<AmbiguousPeptideModSeq, FragmentIon[]> pair : targetPeptides) {
			System.out.println(pair.x.getPeptideAnnotation()+"\t"+pair.y.length);
		}

		// actual localization
		HashSet<FragmentIon> alreadyTaken=new HashSet<FragmentIon>();
		
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=new HashMap<String, Pair<TFloatFloatHashMap,TFloatFloatHashMap>>();
		HashMap<String, HashMap<FragmentIon, XYTrace>> uniqueFragmentIons=new HashMap<String, HashMap<FragmentIon, XYTrace>>();
		HashMap<String, HashMap<FragmentIon, XYTrace>> otherFragmentIons=new HashMap<String, HashMap<FragmentIon, XYTrace>>();
		HashMap<String, FragmentIon[]> uniqueTargetFragments=new HashMap<String, FragmentIon[]>();
		HashMap<String, FragmentIon[]> uniqueIdentifiedTargetFragments=new HashMap<String, FragmentIon[]>();
		HashMap<String, XYPoint> localizationScores=new HashMap<String, XYPoint>();
		
		HashMap<String, TransitionRefinementData> passingForms=new HashMap<String, TransitionRefinementData>();
		TFloatArrayList formsRT=new TFloatArrayList();
		TFloatArrayList scores=new TFloatArrayList(); 
		
		ArrayList<AmbiguousPeptideModSeq> previouslyIdentified=new ArrayList<AmbiguousPeptideModSeq>();
		TFloatArrayList previouslyIdentifiedRTsInSec=new TFloatArrayList(); // can't repeat an RT in sec if an "internal" localization
		FragmentIonBlacklist previouslyIdentifiedIons=new FragmentIonBlacklist(params.getFragmentTolerance());
		
		AmbiguousPeptideModSeq bestPeptideAnnotation=null;
		TransitionRefinementData bestPassingForm=null;
		float bestScore=-Float.MAX_VALUE;
		
		for (Pair<AmbiguousPeptideModSeq, FragmentIon[]> pair : targetPeptides) {
			AmbiguousPeptideModSeq targetPeptideAnnotation=pair.x;
			FragmentIon[] targets=pair.y;
			
			// it's ok that we don't update the targetIons based on previously IDed, for example:
			// if we know we've seen: (S[+80])SSSSK
			// and we're considering: (S[+80]S)SSSK
			// then it's ok if we use matching ions from to (S[+80])SSSSK to identify S(S[+80])SSSK
			// FIXME BUT WE NEED TO REMOVE THE (S[+80])SSSSK RT AS A POSSIBILITY!
			
			// fix ambiguity based on previously identified peptides
			Optional<AmbiguousPeptideModSeq> ambiguityRemoved=targetPeptideAnnotation.removeAmbiguity(previouslyIdentified);
			System.out.println(targetPeptideAnnotation.getPeptideAnnotation()+", "+!ambiguityRemoved.isPresent()); //FIXME
			if (ambiguityRemoved.isPresent()) {
				targetPeptideAnnotation=ambiguityRemoved.get();
			}

			String peptideAnnotation=targetPeptideAnnotation.getPeptideAnnotation();
			String targetPeptideSequence=targetPeptideAnnotation.getPeptideModSeq();
			FragmentationModel model=entryMap.get(targetPeptideSequence);
			if (model==null) {
				// can happen if ambiguity is removed
				model=new FragmentationModel(targetPeptideSequence, params.getAAConstants());
			}
			
			FragmentIon[] allIonsTypes=model.getPrimaryIonObjects(params.getFragType(), precursorCharge);
			double[] allIons=FragmentIon.getMasses(allIonsTypes);
			
			ArrayList<FragmentIon> allTargets=new ArrayList<FragmentIon>(Arrays.asList(targets));
			allTargets.removeAll(alreadyTaken);

			//System.out.println(targetPeptideAnnotation.getPeptideAnnotation()+" --> "+allTargets.size()); //FIXME
			if (allTargets.size()==0) {
				//System.out.println(targetPeptideName+" is degenerate");
				continue;
			}
			targets=allTargets.toArray(new FragmentIon[allTargets.size()]);
			double[] ions=FragmentIon.getMasses(targets);
			//System.out.println("\t"+targetPeptideAnnotation.getPeptideAnnotation()+" --> "+General.toString(ions)); //FIXME
			
			float[] frequencies=background.getFrequencies(ions, precursorMZ, params.getFragmentTolerance());

			//float[] negLogProbsAll=new float[stripes.size()];
			float[] negLogProbsSiteSpecific=new float[stripes.size()];
			for (int k=0; k<stripes.size(); k++) {
				Spectrum spectrum=stripes.get(k);
				//negLogProbsAll[k]=score(params, allIons, allIonsTypes, frequencies, spectrum, false);

				boolean skip=false;
				for (int i=0; i<previouslyIdentifiedRTsInSec.size(); i++) {
					if (spectrum.getScanStartTime()==previouslyIdentifiedRTsInSec.get(i)) {
						// collision!
						skip=true;
						break;
					}
				}
				if (skip) {
					negLogProbsSiteSpecific[k]=-2.0f;
				} else {
					negLogProbsSiteSpecific[k]=score(params, ions, targets, frequencies, spectrum, true);
				}
			}
			//negLogProbsAll=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsAll, Math.round(params.getExpectedPeakWidth()/dutyCycle));//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsAll, movingAverageLength);
			//negLogProbsAll=General.subtract(negLogProbsAll, Log.log10(movingAverageLength)+Log.log10(stripes.size())+Log.log10(peptideModSeqs.size()));
			negLogProbsSiteSpecific=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, Math.round(params.getExpectedPeakWidth()/(dutyCycle)));//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);
			
			negLogProbsSiteSpecific=General.subtract(negLogProbsSiteSpecific, QuickMedian.median(negLogProbsSiteSpecific.clone()));
			
			//negLogProbsSiteSpecific=General.subtract(negLogProbsSiteSpecific, Log.log10(movingAverageLength)+Log.log10(stripes.size())+Log.log10(peptideModSeqs.size()));
			//negLogProbsAll=SkylineSGFilter.paddedSavitzkyGolaySmooth(negLogProbsAll);
			//negLogProbsSiteSpecific=SkylineSGFilter.paddedSavitzkyGolaySmooth(negLogProbsSiteSpecific);

			
			TFloatFloatHashMap coelutingIonsMap=new TFloatFloatHashMap();
			TFloatFloatHashMap uniqueRtScoreMap=new TFloatFloatHashMap();

			if (buildOutGUIData) {
				int[] coelutingIons=TransitionRefiner.numberOfCoelutingIons(ions, allIons, stripes, Math.round((params.getExpectedPeakWidth())/dutyCycle), params.getFragmentTolerance());
				for (int k=0; k<negLogProbsSiteSpecific.length; k++) {
					Spectrum spectrum=stripes.get(k);
					coelutingIonsMap.put(spectrum.getScanStartTime()/60f, coelutingIons[k]);
				}
			}
			
			for (int k=0; k<negLogProbsSiteSpecific.length; k++) {
				Spectrum spectrum=stripes.get(k);
				uniqueRtScoreMap.put(spectrum.getScanStartTime()/60f, negLogProbsSiteSpecific[k]);
			}
			allVsUniqueList.put(peptideAnnotation, new Pair<TFloatFloatHashMap, TFloatFloatHashMap>(coelutingIonsMap, uniqueRtScoreMap));

			EValueCalculator uniqueCalculator=new EValueCalculator(uniqueRtScoreMap);
			float bestRT=uniqueCalculator.getMaxRT()*60f;
			float maxRawScore=uniqueCalculator.getMaxRawScore();
			//System.out.println(targetPeptideAnnotation.getPeptideAnnotation()+" --> "+maxRawScore+" ("+bestRT/60f+")"); //FIXME

			ArrayList<FragmentIon> identifiedTargets=new ArrayList<FragmentIon>();
			Spectrum bestStripe=ChromatogramExtractor.getTargetStripeByRT(stripes, bestRT);
			if (bestStripe!=null) {
				for (int i=0; i<targets.length; i++) {
					float intensity=params.getFragmentTolerance().getIntegratedIntensity(bestStripe.getMassArray(), bestStripe.getIntensityArray(), targets[i].mass);
					if (intensity>0) {
						identifiedTargets.add(targets[i]);
					}
				}
			}

			if (buildOutGUIData) {
				HashMap<FragmentIon, XYTrace> otherTraces=ChromatogramExtractor.extractFragmentChromatograms(params.getFragmentTolerance(), allIonsTypes, stripes, null, GraphType.dashedline);
				HashMap<FragmentIon, XYTrace> uniqueTraces=ChromatogramExtractor.extractFragmentChromatograms(params.getFragmentTolerance(), targets, stripes, null, GraphType.boldline);

				for (FragmentIon ion : uniqueTraces.keySet()) {
					otherTraces.remove(ion);
				}

				uniqueFragmentIons.put(peptideAnnotation, uniqueTraces);
				otherFragmentIons.put(peptideAnnotation, otherTraces);
			}
			uniqueTargetFragments.put(peptideAnnotation, targets);
			uniqueIdentifiedTargetFragments.put(peptideAnnotation, identifiedTargets.toArray(new FragmentIon[identifiedTargets.size()]));

			if (maxRawScore>=minimumScore||maxRawScore>bestScore) {
				bestScore=maxRawScore;
				boolean isLocalized=maxRawScore>=minimumScore&&AmbiguousPeptideModSeq.isLocalized(targetPeptideAnnotation, modification);
				
				if (!AmbiguousPeptideModSeq.isLocalizedAtEnd(targetPeptideAnnotation, modification)) {
					// need to check RTs
					boolean skip=false;
					for (int i=0; i<previouslyIdentifiedRTsInSec.size(); i++) {
						if (bestRT==previouslyIdentifiedRTsInSec.get(i)) {
							// collision!
							skip=true;
							break;
						}
					}
					if (skip) continue;
				}
				
				ArrayList<Spectrum> localStripes=getScanSubset(bestRT-params.getExpectedPeakWidth(), bestRT+params.getExpectedPeakWidth(), stripes);
				TransitionRefinementData quantData=quantifyPeptide(targetPeptideSequence, precursorCharge, targets, bestRT, localStripes, previouslyIdentifiedIons, Optional.ofNullable((float[])null));
				if (quantData==null) continue;
				
				// make sure there are at least 3 "identification" peaks
				int targetNumPeaks=Math.max(params.getMinNumOfQuantitativePeaks(), 3);
				
				int numIdentificationPeaks=0;
				float[] intensities=quantData.getIntegrationArray();
				float[] correlations=quantData.getCorrelationArray();
				FragmentIon[] consideredIons=quantData.getFragmentMassArray();
				ArrayList<FragmentIon> wellShapedIons=new ArrayList<FragmentIon>();
				float localizationIntensity=0.0f;
				for (int i=0; i<correlations.length; i++) {
					if (correlations[i]>=TransitionRefiner.identificationCorrelationThreshold) {
						numIdentificationPeaks++;
						wellShapedIons.add(consideredIons[i]);
						localizationIntensity+=intensities[i];
					}
				}

				if (numIdentificationPeaks==0) {
					// if there's not any localization evidence for a well formed peak then give up
					continue;
				}

				// check for other non-localizing peaks for confirmation
				float[] medianChromatogram=quantData.getMedianChromatogram();
				TransitionRefinementData allQuantData=quantifyPeptide(targetPeptideSequence, precursorCharge, allIonsTypes, bestRT, localStripes, previouslyIdentifiedIons, Optional.of(medianChromatogram));
				if (allQuantData==null) continue;

				float totalIntensity=0.0f;
				numIdentificationPeaks=0;
				intensities=allQuantData.getIntegrationArray();
				correlations=allQuantData.getCorrelationArray();
				for (int i=0; i<correlations.length; i++) {
					if (correlations[i]>=TransitionRefiner.identificationCorrelationThreshold) {
						numIdentificationPeaks++;
						totalIntensity+=intensities[i];
					}
				}
				//System.out.println(targetPeptideAnnotation.getPeptideAnnotation()+" --> "+numIdentificationPeaks+" identification peaks"); //FIXME
				
				// only trust this ID if there are enough peaks!
				if (numIdentificationPeaks>=targetNumPeaks&&quantData.getMedianChromatogram().length>0) {
					
					int numberOfMods=PeptideUtils.getNumberOfMods(targetPeptideSequence, modification.getNominalMass());

					bestRT=quantData.getApexRT();
					formsRT.add(bestRT);

					ModificationLocalizationData modData=new ModificationLocalizationData(targetPeptideAnnotation, bestRT, maxRawScore, numberOfMods, isLocalized, wellShapedIons.toArray(new FragmentIon[wellShapedIons.size()]), localizationIntensity, totalIntensity);

					quantData.setModificationLocalizationData(Optional.of(modData));
					
					bestPassingForm=quantData;
					bestPeptideAnnotation=targetPeptideAnnotation;
					
					if (maxRawScore>minimumScore) {
						alreadyTaken.addAll(Arrays.asList(targets));
						passingForms.put(peptideAnnotation, quantData);
						previouslyIdentified.add(targetPeptideAnnotation);
						
						Range boundaries=allQuantData.getRange();
						float[] rtArray=allQuantData.getRtArray().get();
						for (float f : rtArray) {
							if (boundaries.contains(f)) {
								previouslyIdentifiedRTsInSec.add(f);
							}
						}
						for (FragmentIon target : allIonsTypes) {
							previouslyIdentifiedIons.addIonToBlacklist(target.mass, boundaries);
						}
					}
				}
			}
			localizationScores.put(peptideAnnotation, new XYPoint(bestRT, maxRawScore));

			scores.add(maxRawScore);
		}
		
		if (passingForms.size()==0&&bestPassingForm!=null) {
			passingForms.put(bestPeptideAnnotation.getPeptideAnnotation(), bestPassingForm);
		}
		
		/*if (formsRT.size()==0) {
			System.out.println("multiple\t"+peptideModSeqs.get(0)+"\t"+peptideModSeqs.size()+"\t0\t0\t0");
		} else {
			System.out.println("multiple\t"+peptideModSeqs.get(0)+"\t"+peptideModSeqs.size()+"\t"+formsRT.size()+"\t"+(formsRT.max()-formsRT.min())+"\t"+scores.max());
		}*/
		
		return new PhosphoLocalizationData(allVsUniqueList, uniqueFragmentIons, otherFragmentIons, uniqueTargetFragments, uniqueIdentifiedTargetFragments, localizationScores, passingForms);
	}
	
	public TransitionRefinementData quantifyPeptide(String peptideModSeq, byte precursorCharge, FragmentIon[] targetMasses, float targetRT, ArrayList<Spectrum> stripes, FragmentIonBlacklist blacklistedIons, Optional<float[]> medianChromatogram) {
		/*System.out.println("THINKING ABOUT BLACKLIST: "+peptideModSeq+", "+blacklistedIons.size());
		if (blacklistedIons.size()>0) {
			for (FragmentIon ion : targetMasses) {
				System.out.println("   checking "+ion+" --> "+blacklistedIons.isBlacklisted(ion.mass));
			}
		}*/
		float bestDelta=Float.MAX_VALUE;
		float[] bestIntensities=null;
		ArrayList<float[]> intensityList=new ArrayList<float[]>();
		TFloatArrayList retentionTimes=new TFloatArrayList();
		TFloatArrayList totalIonCurrent=new TFloatArrayList();
		TFloatArrayList totalIdentifiedIonCurrent=new TFloatArrayList();
		for (int k=0; k<stripes.size(); k++) {
			Spectrum spectrum=stripes.get(k);
			retentionTimes.add(spectrum.getScanStartTime());
			float delta=Math.abs(spectrum.getScanStartTime()-targetRT);
			float[] integratedIntensities=params.getFragmentTolerance().getIntegratedIntensities(spectrum.getMassArray(), spectrum.getIntensityArray(), FragmentIon.getMasses(targetMasses));
			for (int i=0; i<integratedIntensities.length; i++) {
				if (integratedIntensities[i]>0.0f&&blacklistedIons.isBlacklisted(targetMasses[i].mass, spectrum.getScanStartTime())) {
						//if (blacklistedIons.size()>0) System.out.println("BLACKLISTED! "+blacklistedIons.isBlacklisted(targetMasses[i].mass, spectrum.getScanStartTime())+"\t"+targetMasses[i]+"\t"+spectrum.getScanStartTime()+"\t"+blacklistedIons.size());
						integratedIntensities[i]=0.0f;
				}
			}
			
			intensityList.add(integratedIntensities);
			if (delta<bestDelta) {
				bestDelta=delta;
				bestIntensities=integratedIntensities;
			}
			float sumIdentifiedIntensities=0.0f;
			for (int i=0; i<integratedIntensities.length; i++) {
				sumIdentifiedIntensities+=integratedIntensities[i];
			}
			totalIdentifiedIonCurrent.add(sumIdentifiedIntensities);
			totalIonCurrent.add(spectrum.getTIC());
		}
		if (bestIntensities.length==0) {
			return null;
		}
		
		// get each scan (fragments by RT)
		TFloatArrayList[] traces=new TFloatArrayList[bestIntensities.length];
		@SuppressWarnings("unchecked")
		ArrayList<XYZPoint>[] deltaMassesByRT=new ArrayList[bestIntensities.length];
		for (int i=0; i<traces.length; i++) {
			traces[i]=new TFloatArrayList();
			deltaMassesByRT[i]=new ArrayList<XYZPoint>();
		}
		for (int index=0; index<intensityList.size(); index++) {
			float[] intensities=intensityList.get(index);
			for (int i=0; i<intensities.length; i++) {
				traces[i].add(intensities[i]);
			}
		}

		// invert each scan into fragment chromatograms (RTs by fragment)
		ArrayList<float[]> chromatograms=new ArrayList<float[]>();
		ArrayList<ArrayList<XYZPoint>> chromatogramDeltaMassesByRT=new ArrayList<ArrayList<XYZPoint>>();
		TFloatArrayList keptIntensities=new TFloatArrayList();
		ArrayList<FragmentIon> keptMasses=new ArrayList<FragmentIon>();
		for (int i=0; i<bestIntensities.length; i++) {
			if (bestIntensities[i]>0.0f) {
				float[] chromatogram=traces[i].toArray();
				chromatogram=SkylineSGFilter.paddedSavitzkyGolaySmooth(chromatogram);
				chromatograms.add(chromatogram);
				chromatogramDeltaMassesByRT.add(deltaMassesByRT[i]);
				keptIntensities.add(bestIntensities[i]);
				keptMasses.add(targetMasses[i]);
			}
		}

		// identify transitions
		TransitionRefinementData data=TransitionRefiner.identifyTransitions(peptideModSeq, precursorCharge, keptMasses.toArray(new FragmentIon[keptMasses.size()]), chromatograms, retentionTimes.toArray(), medianChromatogram);
		float[] correlations=data.getCorrelationArray();
		float[] integrations=data.getIntegrationArray();
		Range rtRange=data.getRange();

		TDoubleArrayList mzs=new TDoubleArrayList();
		TFloatArrayList intens=new TFloatArrayList();
		TFloatArrayList deltaMasses=new TFloatArrayList(); // will ultimately be the length of the correlations array

		float correlationThreshold=-1f; // Keep all localizing fragments!
		for (int i=0; i<keptIntensities.size(); i++) {
			// calculate delta mass for each fragment ion
			
			float totalDeltaMasses=0.0f;
			float totalIntensities=0.0f;
			for (XYZPoint point : chromatogramDeltaMassesByRT.get(i)) {
				if (rtRange.contains((float)point.getX())) {
					totalDeltaMasses+=point.getY();
					totalIntensities+=point.getZ();
				}
			}
			
			float deltaMass=0.0f;
			if (totalIntensities>0.0f) {
				deltaMass=totalDeltaMasses/totalIntensities;
			}
			deltaMasses.add(deltaMass);
			
			// generate spectrum peaks for only each "kept" fragment ion
			if (correlations[i]>=correlationThreshold) {
				// grab mz
				if (keptIntensities.get(i)>0) {
					mzs.add(keptMasses.get(i).mass);
					intens.add(integrations[i]);
				}
			}
		}

		float ticSum=0.0f;
		float identifiedTicSum=0.0f;
		for (int i=0; i<stripes.size(); i++) {
			Spectrum stripe=stripes.get(i);
			if (rtRange.contains(stripe.getScanStartTime())) {
				identifiedTicSum+=totalIdentifiedIonCurrent.get(i);
				ticSum+=totalIonCurrent.get(i);
			}
		}
		float identifiedTICRatio=ticSum==0.0f?0.0f:identifiedTicSum/ticSum;

		double[] massArray=mzs.toArray();
		float[] intensityArray=intens.toArray();
		float[] deltaMassArray=deltaMasses.toArray();
		return data.addPeakData(deltaMassArray, massArray, intensityArray, retentionTimes.toArray(), identifiedTICRatio);
	}
	
	public static float score(SearchParameters parameters, double[] ions, FragmentIon[] ionTypes, float[] frequencies, Spectrum stripe, boolean report) {
		if (frequencies.length==0) return 0.0f;

		double[] massArray=stripe.getMassArray();
		
		// only keep 
		HashSet<String> uniqueIonTypes=new HashSet<String>();
		TObjectFloatHashMap<String> matches=new TObjectFloatHashMap<String>();
		
		for (int i=0; i<frequencies.length; i++) {
			String canonicalIonType=ionTypes[i].toCanonicalIonTypeString();
			uniqueIonTypes.add(canonicalIonType);
			
			boolean match=parameters.getFragmentTolerance().getIndex(massArray, ions[i]).isPresent();
			if (match) {
				if (matches.contains(canonicalIonType)) {
					if (matches.get(canonicalIonType)>frequencies[i]) {
						matches.put(canonicalIonType, frequencies[i]);
					}
				} else {
					matches.put(canonicalIonType, frequencies[i]);
				}
			}
		}

		float logProb=Log.log10(uniqueIonTypes.size()); // bonferroni correction
		for (float freq : matches.values()) {
			logProb+=Log.protectedLog10(freq);
		}

		return -logProb;
	}

	public static FragmentIon[] getUniqueFragmentIons(String peptideModSeq, byte precursorCharge, HashMap<String, FragmentationModel> availableModels, SearchParameters params) {
		FragmentationModel unitEntry=availableModels.get(peptideModSeq);
		HashSet<FragmentIon> ions=new HashSet<FragmentIon>(Arrays.asList(unitEntry.getPrimaryIonObjects(params.getFragType(), precursorCharge, false)));

		for (Entry<String, FragmentationModel> otherEntry : availableModels.entrySet()) {
			String otherPeptideModSeq=otherEntry.getKey();
			if (!peptideModSeq.equals(otherPeptideModSeq)) {
				FragmentationModel otherUnitEntry=otherEntry.getValue();
				ions.removeAll(Arrays.asList(otherUnitEntry.getPrimaryIonObjects(params.getFragType(), precursorCharge, false)));
			}
		}

		FragmentIon[] ionArray=ions.toArray(new FragmentIon[ions.size()]);
		Arrays.sort(ionArray);
		return ionArray;
	}

	public static HashMap<String, FragmentIon[]> getUniqueFragmentIons(byte precursorCharge, HashMap<String, FragmentationModel> entryMap, SearchParameters params) {
		HashMap<String, FragmentIon[]> uniqueIons=new HashMap<String, FragmentIon[]>();
		for (Entry<String, FragmentationModel> entry : entryMap.entrySet()) {
			String peptideModSeq=entry.getKey();
			FragmentationModel unitEntry=entry.getValue();
			HashSet<FragmentIon> ions=new HashSet<FragmentIon>(Arrays.asList(unitEntry.getPrimaryIonObjects(params.getFragType(), precursorCharge, false)));

			for (Entry<String, FragmentationModel> otherEntry : entryMap.entrySet()) {
				String otherPeptideModSeq=otherEntry.getKey();
				if (peptideModSeq!=otherPeptideModSeq) {
					// actual != is ok here because we're dealing with the same objects
					FragmentationModel otherUnitEntry=otherEntry.getValue();
					ions.removeAll(Arrays.asList(otherUnitEntry.getPrimaryIonObjects(params.getFragType(), precursorCharge, false)));
				}
			}
			
			if (ions.size()>0) {
				FragmentIon[] ionArray=ions.toArray(new FragmentIon[ions.size()]);
				Arrays.sort(ionArray);
				uniqueIons.put(peptideModSeq, ionArray);
			}
		}
		return uniqueIons;
	}
	
	public static ArrayList<Spectrum> getScanSubsetFromStripes(float minRT, float maxRT, ArrayList<Stripe> allScansInStripe) {
		ArrayList<Spectrum> subset=new ArrayList<Spectrum>();
		for (Spectrum scan : allScansInStripe) {
			if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
				subset.add(scan);
			}
		}
		return subset;
	}
	
	public static ArrayList<Spectrum> getScanSubset(float minRT, float maxRT, ArrayList<Spectrum> allScansInStripe) {
		ArrayList<Spectrum> subset=new ArrayList<Spectrum>();
		for (Spectrum scan : allScansInStripe) {
			if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
				subset.add(scan);
			}
		}
		return subset;
	}
}
