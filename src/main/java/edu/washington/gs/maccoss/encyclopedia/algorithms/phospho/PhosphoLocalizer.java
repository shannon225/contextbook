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
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
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
			return new PhosphoLocalizationData(new HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>>(), new HashMap<String, XYTrace[]>(), new HashMap<String, XYPoint>());
		} else {
			PhosphoLocalizationData multiple=extractPhosphoForms(psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), permutations, psmdata.getRetentionTime(), stripes);
			return multiple;
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
	PhosphoLocalizationData extractPhosphoForms(double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, float retentionTime, ArrayList<Spectrum> allScansInStripe) {
		float dutyCycle=1.0f;
		for (Entry<Range, Float> entry : diaFile.getRanges().entrySet()) {
			if (entry.getKey().contains((float)precursorMZ)) {
				dutyCycle=entry.getValue();
				break;
			}
		}
		int movingAverageLength=Math.round(params.getExpectedPeakWidth()/dutyCycle/2.0f);
		
		float duration=5f*60f; // search for 5 minutes

		ArrayList<Spectrum> stripes=getScanSubset(retentionTime-duration, retentionTime+duration, allScansInStripe);
		
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
		
		ArrayList<Triplet<String, String, FragmentIon[]>> targetPeptidesLeft=new ArrayList<Triplet<String,String,FragmentIon[]>>(); 
		// go left to right, drop the last
		for (int i=0; i<peptideModSeqs.size()-1; i++) {
			String targetPeptide=peptideModSeqs.get(i);
			String targetPeptideName=getLeftAnnotation(targetPeptide);

			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of
			// remaining sequences to the right
			for (int j=peptideModSeqs.size()-1; j>=i; j--) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			FragmentIon[] targets=PhosphoLocalizer.getUniqueFragmentIons(targetPeptide, precursorCharge, modelBatch, params);
			targetPeptidesLeft.add(new Triplet<String, String, FragmentIon[]>(targetPeptide, targetPeptideName, targets));
		}
		ArrayList<Triplet<String, String, FragmentIon[]>> targetPeptidesRight=new ArrayList<Triplet<String,String,FragmentIon[]>>(); 
		// go right to left, drop the first
		for (int i=peptideModSeqs.size()-1; i>=1; i--) {
			String targetPeptide=peptideModSeqs.get(i);
			String targetPeptideName=getRightAnnotation(targetPeptide);
			
			HashMap<String, FragmentationModel> modelBatch=new HashMap<String, FragmentationModel>();
			// shrink the number of unique ions subtractors to the pool of remaining sequences to the right
			for (int j=peptideModSeqs.size()-1; j>=i; j--) {
				String seq=peptideModSeqs.get(j);
				modelBatch.put(seq, entryMap.get(seq));
			}
			FragmentIon[] targets=PhosphoLocalizer.getUniqueFragmentIons(targetPeptide, precursorCharge, modelBatch, params);
			targetPeptidesRight.add(new Triplet<String, String, FragmentIon[]>(targetPeptide, targetPeptideName, targets));
		}
		
		// interlay the peptides so we look at the most localizing first
		ArrayList<Triplet<String, String, FragmentIon[]>> targetPeptides=new ArrayList<Triplet<String,String,FragmentIon[]>>();
		for (int i=0; i<targetPeptidesLeft.size(); i++) {
			targetPeptides.add(targetPeptidesLeft.get(i));
			targetPeptides.add(targetPeptidesRight.get(i));
		}

		// actual localization
		HashSet<FragmentIon> alreadyTaken=new HashSet<FragmentIon>();
		
		HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> allVsUniqueList=new HashMap<String, Pair<TFloatFloatHashMap,TFloatFloatHashMap>>();
		HashMap<String, XYTrace[]> uniqueFragmentIons=new HashMap<String, XYTrace[]>();
		HashMap<String, XYPoint> localizationScores=new HashMap<String, XYPoint>();
		
		TFloatArrayList formsRT=new TFloatArrayList();
		TFloatArrayList scores=new TFloatArrayList(); 
		
		for (Triplet<String, String, FragmentIon[]> triplet : targetPeptides) {
			String targetPeptide=triplet.x;
			String targetPeptideName=triplet.y;
			FragmentIon[] targets=triplet.z;
			
			FragmentationModel model=entryMap.get(targetPeptide);
			FragmentIon[] allIonsTypes=model.getPrimaryIonObjects(params.getFragType(), precursorCharge);
			double[] allIons=FragmentIon.getMasses(allIonsTypes);
			
			ArrayList<FragmentIon> allTargets=new ArrayList<FragmentIon>(Arrays.asList(targets));
			allTargets.removeAll(alreadyTaken);
			if (allTargets.size()==0) {
				//System.out.println(targetPeptideName+" is degenerate");
				continue;
			}
			targets=allTargets.toArray(new FragmentIon[allTargets.size()]);
			double[] ions=FragmentIon.getMasses(targets);
			
			float[] frequencies=background.getFrequencies(ions, precursorMZ, params.getFragmentTolerance());

			float[] negLogProbsAll=new float[stripes.size()];
			float[] negLogProbsSiteSpecific=new float[stripes.size()];
			for (int k=0; k<stripes.size(); k++) {
				Spectrum spectrum=stripes.get(k);
				negLogProbsAll[k]=score(params, allIons, allIonsTypes, frequencies, spectrum, false);
				negLogProbsSiteSpecific[k]=score(params, ions, targets, frequencies, spectrum, true);
			}
			negLogProbsAll=AbstractLibraryScoringTask.movingCenteredSum(negLogProbsAll, movingAverageLength);//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsAll, movingAverageLength);
			negLogProbsAll=General.subtract(negLogProbsAll, Log.log10(movingAverageLength)+Log.log10(stripes.size())+Log.log10(peptideModSeqs.size()));
			negLogProbsSiteSpecific=AbstractLibraryScoringTask.movingCenteredSum(negLogProbsSiteSpecific, movingAverageLength);//AbstractLibraryScoringTask.gaussianCenteredAverage(negLogProbsSiteSpecific, movingAverageLength);
			negLogProbsSiteSpecific=General.subtract(negLogProbsSiteSpecific, Log.log10(movingAverageLength)+Log.log10(stripes.size())+Log.log10(peptideModSeqs.size()));

			TFloatFloatHashMap allRtScoreMap=new TFloatFloatHashMap();
			TFloatFloatHashMap uniqueRtScoreMap=new TFloatFloatHashMap();
			for (int k=0; k<negLogProbsSiteSpecific.length; k++) {
				Spectrum spectrum=stripes.get(k);
				allRtScoreMap.put(spectrum.getScanStartTime()/60f, negLogProbsAll[k]);
				uniqueRtScoreMap.put(spectrum.getScanStartTime()/60f, negLogProbsSiteSpecific[k]);
			}
			allVsUniqueList.put(targetPeptideName, new Pair<TFloatFloatHashMap, TFloatFloatHashMap>(allRtScoreMap, uniqueRtScoreMap));

			EValueCalculator uniqueCalculator=new EValueCalculator(uniqueRtScoreMap);

			XYTrace[] traces=ChromatogramExtractor.extractFragmentChromatograms(params.getFragmentTolerance(), targets, stripes);
			uniqueFragmentIons.put(targetPeptideName, traces);
			//Charter.launchChart("Retention Time (Site Specific)", "Intensity", false, new Dimension(800, 250), traces);
			//traces=ChromatogramExtractor.extractFragmentChromatograms(params.getFragmentTolerance(), totalIons.toArray(new FragmentIon[totalIons.size()]), stripes);
			//Charter.launchChart("Retention Time (All Ions)", "Intensity", false, new Dimension(800, 250), traces);
			
			//Charter.launchChart("All Score", "Count", true, allCalculator.toTraces());
			//Charter.launchChart("Unique Score", "Count", true, uniqueCalculator.toTraces());

			float bestRT=uniqueCalculator.getMaxRT();
			float maxRawScore=uniqueCalculator.getMaxRawScore();
			localizationScores.put(targetPeptideName, new XYPoint(bestRT, maxRawScore));
			if (maxRawScore>2f) {
				formsRT.add(bestRT);
				alreadyTaken.addAll(Arrays.asList(targets));
				//System.out.println(targetPeptideName+" kept, score: "+uniqueCalculator.getMaxRawScore());
			}

			scores.add(maxRawScore);
		}
		
		if (formsRT.size()==0) {
			System.out.println("multiple\t"+peptideModSeqs.get(0)+"\t"+peptideModSeqs.size()+"\t0\t0\t0");
		} else {
			System.out.println("multiple\t"+peptideModSeqs.get(0)+"\t"+peptideModSeqs.size()+"\t"+formsRT.size()+"\t"+(formsRT.max()-formsRT.min())+"\t"+scores.max());
		}
		
		return new PhosphoLocalizationData(allVsUniqueList, uniqueFragmentIons, localizationScores);
	}

	public static String getLeftAnnotation(String targetPeptide) {
		StringBuilder sb=new StringBuilder(targetPeptide);
		for (int j=0; j<sb.length(); j++) {
			char c=sb.charAt(j);
			if (c=='S'||c=='T'||c=='Y') {
				sb.insert(j, "(");
				break;
			}
		}
		
		String phospho="[+79.96633]";
		int index=sb.lastIndexOf(phospho);
		sb.insert(index+phospho.length(), ")");
		String targetPeptideName=sb.toString();
		return targetPeptideName;
	}

	public static String getRightAnnotation(String targetPeptide) {
		StringBuilder sb=new StringBuilder(targetPeptide);
		String phospho="[+79.96633]";
		int index=sb.indexOf(phospho);
		sb.insert(index-1, "(");

		int lastIndex=sb.lastIndexOf(phospho);
		for (int j=sb.length()-1; j>=0; j--) {
			if (j<lastIndex+phospho.length()) {
				sb.insert(lastIndex+phospho.length(), ")");
				break;
			}
			char c=sb.charAt(j);
			if (c=='S'||c=='T'||c=='Y') {
				sb.insert(j+1, ")");
				break;
			}
		}
		String targetPeptideName=sb.toString();
		return targetPeptideName;
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

	public static FragmentIon[] getUniqueFragmentIons(String peptideModSeq, byte precursorCharge, HashMap<String, FragmentationModel> availableModels, SearchParameters params) {
		FragmentationModel unitEntry=availableModels.get(peptideModSeq);
		HashSet<FragmentIon> ions=new HashSet<FragmentIon>(Arrays.asList(unitEntry.getPrimaryIonObjects(params.getFragType(), precursorCharge)));

		for (Entry<String, FragmentationModel> otherEntry : availableModels.entrySet()) {
			String otherPeptideModSeq=otherEntry.getKey();
			if (!peptideModSeq.equals(otherPeptideModSeq)) {
				FragmentationModel otherUnitEntry=otherEntry.getValue();
				ions.removeAll(Arrays.asList(otherUnitEntry.getPrimaryIonObjects(params.getFragType(), precursorCharge)));
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
