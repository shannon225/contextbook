package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

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
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
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

	public HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> runPhosphoLocalization(PSMData psmdata, ArrayList<Stripe> stripes) {
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(psmdata.getPeptideModSeq(), params.getAAConstants());
		if (permutations.size()==1) {
			System.out.println("single\t"+psmdata.getPeptideModSeq());
			return new HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>>();
		} else {
			HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> multiple=extractPhosphoForms(psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), permutations, psmdata.getRetentionTime(), stripes);
			System.out.println("multiple\t"+psmdata.getPeptideModSeq()+"\t"+multiple);
			return multiple;
		}
	}

	private HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> extractPhosphoForms(double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, float retentionTime, ArrayList<Stripe> allScansInStripe) {
		float dutyCycle=1.0f;
		for (Entry<Range, Float> entry : diaFile.getRanges().entrySet()) {
			if (entry.getKey().contains((float)precursorMZ)) {
				dutyCycle=entry.getValue();
				break;
			}
		}
		int movingAverageLength=Math.round(params.getExpectedPeakWidth()/dutyCycle);
		
		float duration=6*60f; // search for 6 minutes

		ArrayList<Stripe> stripes=allScansInStripe;//getScanSubset(retentionTime-duration, retentionTime+duration, allScansInStripe);
		
		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
		for (String peptideModSeq : peptideModSeqs) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			entryMap.put(peptideModSeq, model);
		}
		
		HashMap<String, FragmentIon[]> uniqueIons=getUniqueFragmentIons(precursorCharge, entryMap, params);
		
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=new HashMap<String, Pair<TFloatFloatHashMap,TFloatFloatHashMap>>();
		for (Entry<String, FragmentationModel> entry : entryMap.entrySet()) {
			String peptideModSeq=entry.getKey();
			FragmentationModel model=entry.getValue();
			double[] allIons=model.getPrimaryIons(params.getFragType(), precursorCharge);
			
			FragmentIon[] targets=uniqueIons.get(peptideModSeq);
			double[] ions=FragmentIon.getMasses(targets);
			float[] frequencies=background.getFrequencies(ions, precursorMZ, params.getFragmentTolerance());

			float[] negLogProbsAll=new float[stripes.size()];
			float[] negLogProbsSiteSpecific=new float[stripes.size()];
			for (int i=0; i<stripes.size(); i++) {
				Stripe spectrum=stripes.get(i);
				negLogProbsAll[i]=score(params, allIons, frequencies, spectrum);
				negLogProbsSiteSpecific[i]=score(params, ions, frequencies, spectrum);
			}
			negLogProbsAll=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsAll, movingAverageLength);
			negLogProbsSiteSpecific=AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);

			TFloatFloatHashMap allRtScoreMap=new TFloatFloatHashMap();
			TFloatFloatHashMap uniqueRtScoreMap=new TFloatFloatHashMap();
			for (int i=0; i<negLogProbsSiteSpecific.length; i++) {
				Stripe spectrum=stripes.get(i);
				allRtScoreMap.put(spectrum.getScanStartTime()/60f, negLogProbsAll[i]);
				uniqueRtScoreMap.put(spectrum.getScanStartTime()/60f, negLogProbsSiteSpecific[i]);
			}
			allVsUniqueList.put(peptideModSeq, new Pair<TFloatFloatHashMap, TFloatFloatHashMap>(allRtScoreMap, uniqueRtScoreMap));

			EValueCalculator allCalculator=new EValueCalculator(allRtScoreMap);
			EValueCalculator uniqueCalculator=new EValueCalculator(uniqueRtScoreMap);

			float bestRT=uniqueCalculator.getMaxRT();
			float allScore=allRtScoreMap.get(bestRT);
			System.out.println("FINAL: "+peptideModSeq+" --> "+bestRT/60.0f+"/"+allCalculator.getMaxRT()/60.0f+", site specific: "+uniqueCalculator.getMaxRawScore()+" ("+uniqueCalculator.getNegLog10EValue(bestRT)+"), all: "+allScore+" ("+allCalculator.getNegLog10EValue(allScore)+")");
		}
		return allVsUniqueList;
	}
	
	private static float score(SearchParameters parameters, double[] ions, float[] frequencies, Stripe stripe) {
		if (frequencies.length==0) return 0.0f;

		double[] massArray=stripe.getMassArray();
		float logProb=0.0f;
		for (int i=0; i<frequencies.length; i++) {
			float hitProb=frequencies[i]*massArray.length;
			boolean match=parameters.getFragmentTolerance().getIndex(massArray, ions[i]).isPresent();
			if (match) {
				logProb+=Log.log10(hitProb);
			}
		}
		// neg log prob (normalized by N attempts)
		return Log.log10(frequencies.length)-logProb;
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
			FragmentIon[] ionArray=ions.toArray(new FragmentIon[ions.size()]);
			Arrays.sort(ionArray);
			uniqueIons.put(peptideModSeq, ionArray);
		}
		return uniqueIons;
	}
	
	public ArrayList<Stripe> getScanSubset(float minRT, float maxRT, ArrayList<Stripe> allScansInStripe) {
		ArrayList<Stripe> subset=new ArrayList<Stripe>();
		for (Stripe scan : allScansInStripe) {
			if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
				subset.add(scan);
			}
		}
		return subset;
	}
}
