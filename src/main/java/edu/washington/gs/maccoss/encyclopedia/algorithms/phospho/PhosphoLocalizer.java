package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.awt.Dimension;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.ChromatogramExtractor;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PhosphoLocalizer {
	private final StripeFileInterface diaFile;
	private final SearchParameters params;
	private final BackgroundFrequencyCalculator background;

	public PhosphoLocalizer(StripeFileInterface diaFile, LibraryInterface searchedLibrary, SearchParameters params) throws IOException,DataFormatException,SQLException {
		this.diaFile=diaFile;
		this.params=params;
		background=BackgroundFrequencyCalculator.generateBackground(diaFile, searchedLibrary);
	}

	public PhosphoLocalizationData runDIAPhosphoLocalization(PSMData psmdata, ArrayList<Stripe> stripes) {
		ArrayList<Spectrum> spectra=new ArrayList<Spectrum>();
		for (Stripe stripe : stripes) {
			spectra.add(stripe);
		}
		return runPhosphoLocalization(psmdata, spectra);
	}

	public PhosphoLocalizationData runPhosphoLocalization(PSMData psmdata, ArrayList<Spectrum> stripes) {
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(psmdata.getPeptideModSeq(), params.getAAConstants());
		if (permutations.size()==1) {
			System.out.println("single\t"+psmdata.getPeptideModSeq()+"\t1\t1\t0\t1000");
			return new PhosphoLocalizationData(new HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>>());
		} else {
			PhosphoLocalizationData multiple=extractPhosphoForms(psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), permutations, psmdata.getRetentionTime(), stripes);
			return multiple;
		}
	}

	PhosphoLocalizationData extractPhosphoForms(double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, float retentionTime, ArrayList<Spectrum> allScansInStripe) {
		float dutyCycle=1.0f;
		for (Entry<Range, Float> entry : diaFile.getRanges().entrySet()) {
			if (entry.getKey().contains((float)precursorMZ)) {
				dutyCycle=entry.getValue();
				break;
			}
		}
		int movingAverageLength=Math.round(params.getExpectedPeakWidth()/dutyCycle/2.0f);
		
		float duration=1.5f*60f; // search for 5 minutes

		ArrayList<Spectrum> stripes=getScanSubset(retentionTime-duration, retentionTime+duration, allScansInStripe);
		
		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
		for (String peptideModSeq : peptideModSeqs) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			entryMap.put(peptideModSeq, model);
		}

		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=new HashMap<String, Pair<TFloatFloatHashMap,TFloatFloatHashMap>>();
		HashSet<FragmentIon> alreadyTaken=new HashSet<FragmentIon>();
		TFloatArrayList formsRT=new TFloatArrayList();
		TFloatArrayList scores=new TFloatArrayList(); 
		scores.add(0.0f);
		int round=0;
		while (!entryMap.isEmpty()) {
			round++;
			if (round>1) break;
			
			HashMap<String, FragmentIon[]> uniqueIons=getUniqueFragmentIons(precursorCharge, entryMap, params);
			if (uniqueIons.size()==0) {
				break;
			}
			
			if (uniqueIons.size()==1&&alreadyTaken.size()==0) {
				// can't discriminate between forms at all
				//break;
			}
			HashSet<FragmentIon> totalIons=new HashSet<FragmentIon>();

			for (FragmentationModel model : entryMap.values()) {
				FragmentIon[] allIonsTypes=model.getPrimaryIonObjects(params.getFragType(), precursorCharge);
				totalIons.addAll(Arrays.asList(allIonsTypes));
			}

			for (Entry<String, FragmentIon[]> entry : uniqueIons.entrySet()) {
				String peptideModSeq=entry.getKey();
				FragmentIon[] targets=entry.getValue();
				FragmentationModel model=entryMap.remove(peptideModSeq);
				FragmentIon[] allIonsTypes=model.getPrimaryIonObjects(params.getFragType(), precursorCharge);
				double[] allIons=FragmentIon.getMasses(allIonsTypes);
				
				ArrayList<FragmentIon> allTargets=new ArrayList<FragmentIon>(Arrays.asList(targets));
				allTargets.removeAll(alreadyTaken);
				if (allTargets.size()==0) {
					continue;
				}
				targets=allTargets.toArray(new FragmentIon[allTargets.size()]);
				double[] ions=FragmentIon.getMasses(targets);

				/*System.out.println(peptideModSeq+" +"+precursorCharge);
				for (FragmentIon ion : targets) {
					System.out.println("\t"+ion+", "+Math.round(ion.mass));
				}*/
				
				float[] frequencies=background.getFrequencies(ions, precursorMZ, params.getFragmentTolerance());

				float[] negLogProbsAll=new float[stripes.size()];
				float[] negLogProbsSiteSpecific=new float[stripes.size()];
				for (int i=0; i<stripes.size(); i++) {
					Spectrum spectrum=stripes.get(i);
					negLogProbsAll[i]=score(params, allIons, allIonsTypes, frequencies, spectrum, false);
					negLogProbsSiteSpecific[i]=score(params, ions, targets, frequencies, spectrum, true);
				}
				negLogProbsAll=AbstractLibraryScoringTask.movingCenteredSum(negLogProbsAll, movingAverageLength);//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsAll, movingAverageLength);
				negLogProbsAll=General.subtract(negLogProbsAll, Log.log10(movingAverageLength)+Log.log10(peptideModSeqs.size()));
				negLogProbsSiteSpecific=AbstractLibraryScoringTask.movingCenteredSum(negLogProbsSiteSpecific, movingAverageLength);//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);
				negLogProbsSiteSpecific=General.subtract(negLogProbsSiteSpecific, Log.log10(movingAverageLength)+Log.log10(peptideModSeqs.size()));

				TFloatFloatHashMap allRtScoreMap=new TFloatFloatHashMap();
				TFloatFloatHashMap uniqueRtScoreMap=new TFloatFloatHashMap();
				for (int i=0; i<negLogProbsSiteSpecific.length; i++) {
					Spectrum spectrum=stripes.get(i);
					allRtScoreMap.put(spectrum.getScanStartTime()/60f, negLogProbsAll[i]);
					uniqueRtScoreMap.put(spectrum.getScanStartTime()/60f, negLogProbsSiteSpecific[i]);
				}
				allVsUniqueList.put(peptideModSeq, new Pair<TFloatFloatHashMap, TFloatFloatHashMap>(allRtScoreMap, uniqueRtScoreMap));

				EValueCalculator uniqueCalculator=new EValueCalculator(uniqueRtScoreMap);

				XYTrace[] traces=ChromatogramExtractor.extractFragmentChromatograms(params.getFragmentTolerance(), targets, stripes);
				Charter.launchChart("Retention Time (Site Specific)", "Intensity", false, new Dimension(800, 250), traces);
				traces=ChromatogramExtractor.extractFragmentChromatograms(params.getFragmentTolerance(), totalIons.toArray(new FragmentIon[totalIons.size()]), stripes);
				Charter.launchChart("Retention Time (All Ions)", "Intensity", false, new Dimension(800, 250), traces);
				
				//Charter.launchChart("All Score", "Count", true, allCalculator.toTraces());
				//Charter.launchChart("Unique Score", "Count", true, uniqueCalculator.toTraces());

				float bestRT=uniqueCalculator.getMaxRT();
				if (uniqueCalculator.getMaxRawScore()>2f) {
					formsRT.add(bestRT);
					alreadyTaken.addAll(Arrays.asList(targets));
				}
				if (round==1) {
					scores.add(uniqueCalculator.getMaxRawScore());
				}
				//System.out.println("Score: "+uniqueCalculator.getMaxRawScore()+"\n");
			}
		}
		if (formsRT.size()==0) {
			System.out.println("multiple\t"+peptideModSeqs.get(0)+"\t"+peptideModSeqs.size()+"\t0\t0\t0");
		} else {
			System.out.println("multiple\t"+peptideModSeqs.get(0)+"\t"+peptideModSeqs.size()+"\t"+formsRT.size()+"\t"+(formsRT.max()-formsRT.min())+"\t"+scores.max());
		}
		
		return new PhosphoLocalizationData(allVsUniqueList);
	}
	
	private static float score(SearchParameters parameters, double[] ions, FragmentIon[] ionTypes, float[] frequencies, Spectrum stripe, boolean report) {
		if (frequencies.length==0) return 0.0f;

		double[] massArray=stripe.getMassArray();
		float logProb=0.0f;
		ArrayList<FragmentIon> matches=new ArrayList<FragmentIon>();
		for (int i=0; i<frequencies.length; i++) {
			float hitProb=frequencies[i]*massArray.length;
			boolean match=parameters.getFragmentTolerance().getIndex(massArray, ions[i]).isPresent();
			if (match) {
				logProb+=Log.log10(hitProb);
				matches.add(ionTypes[i]);
			}
		}
		//if (report&&matches.size()>0) System.out.println(stripe.getScanStartTime()/60f+"\tFound:"+General.toString(matches)+" ("+matches.size()+"/"+ions.length+")\t"+(-logProb-Log.log10(frequencies.length)));
		// neg log prob (normalized by N attempts)
		return -logProb-Log.log10(frequencies.length);
	}

	public static HashMap<String, FragmentIon[]> getUniqueFragmentIons(byte precursorCharge, HashMap<String, FragmentationModel> entryMap, SearchParameters params) {
		HashMap<String, FragmentIon[]> uniqueIons=new HashMap<String, FragmentIon[]>();
		for (Entry<String, FragmentationModel> entry : entryMap.entrySet()) {
			String peptideModSeq=entry.getKey();
			FragmentationModel unitEntry=entry.getValue();
			HashSet<FragmentIon> ions=new HashSet<FragmentIon>(Arrays.asList(unitEntry.getPrimaryIonObjects(params.getFragType(), precursorCharge)));

			for (Entry<String, FragmentationModel> otherEntry : entryMap.entrySet()) {
				String otherPeptideModSeq=otherEntry.getKey();
				if (peptideModSeq!=otherPeptideModSeq) {
					// actual != is ok here because we're dealing with the same objects
					FragmentationModel otherUnitEntry=otherEntry.getValue();
					ions.removeAll(Arrays.asList(otherUnitEntry.getPrimaryIonObjects(params.getFragType(), precursorCharge)));
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
	
	public ArrayList<Spectrum> getScanSubset(float minRT, float maxRT, ArrayList<Spectrum> allScansInStripe) {
		ArrayList<Spectrum> subset=new ArrayList<Spectrum>();
		for (Spectrum scan : allScansInStripe) {
			if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
				subset.add(scan);
			}
		}
		return subset;
	}
}
