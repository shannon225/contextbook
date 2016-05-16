package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakScores;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.map.hash.TFloatObjectHashMap;
import gnu.trove.procedure.TFloatObjectProcedure;
import gnu.trove.set.hash.TDoubleHashSet;
import gnu.trove.set.hash.TFloatHashSet;

public class PhosphoLocalizer {
	private static final String unknownKey="unknown";
	private final String filename;
	private final SearchParameters params;
	private final PSMData psmdata;
	private final ArrayList<Stripe> stripes;

	public PhosphoLocalizer(String filename, SearchParameters params, PSMData psmdata, ArrayList<Stripe> stripes) {
		this.filename=filename;
		this.params=params;
		this.psmdata=psmdata;
		this.stripes=stripes;
	}

	public void runPhosphoLocalization() {
		ArrayList<String> permutations=PhosphoPermuter.getPermutations(psmdata.getPeptideModSeq(), params.getAAConstants());
		if (permutations.size()==1) {
			System.out.println("single\t"+psmdata.getPeptideModSeq());
		} else {
			boolean multiple=extractPhosphoForms(psmdata.getPrecursorMZ(), psmdata.getPrecursorCharge(), permutations, psmdata.getRetentionTime());
			System.out.println("multiple\t"+psmdata.getPeptideModSeq()+"\t"+multiple);
		}
	}

	private boolean extractPhosphoForms(double precursorMZ, byte precursorCharge, ArrayList<String> peptideModSeqs, float retentionTime) {
		float duration=6*60f; // search for 6 minutes
		
		EncyclopediaOneScorer encyclopediaScorer=new EncyclopediaOneScorer(params, null); // not using aux scoring

		ArrayList<Stripe> stripes=getScanSubset(retentionTime-duration, retentionTime+duration);
		
		HashMap<String, FragmentationModel> entryMap=new HashMap<String, FragmentationModel>();
		for (String peptideModSeq : peptideModSeqs) {
			FragmentationModel model=new FragmentationModel(peptideModSeq, params.getAAConstants());
			entryMap.put(peptideModSeq, model);
		}
		
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

		final TFloatObjectHashMap<String> allBestTimes=new TFloatObjectHashMap<String>();
		TFloatHashSet list=new TFloatHashSet();
		
		for (Entry<String, FragmentationModel> entry : entryMap.entrySet()) {
			String peptideModSeq=entry.getKey();
			LibraryEntry unitEntry=entry.getValue().getUnitSpectrum(filename, new HashSet<String>(), precursorCharge, retentionTime, params);

			TFloatObjectHashMap<String> bestTimes=new TFloatObjectHashMap<String>();
			float bestScore=0.0f;
			TFloatFloatHashMap rtScoreMap=new TFloatFloatHashMap();
			for (Stripe spectrum : stripes) {
				// float score=encyclopediaScorer.score(unitEntry, spectrum);
				PeakScores[] individualPeakScores=encyclopediaScorer.getIndividualPeakScores(unitEntry, spectrum, false);
				float score=score(individualPeakScores);
				
				rtScoreMap.put(spectrum.getScanStartTime(), score);
				if (score>9.5f) { // at least 5 transitions
					PeakScores[] uniquePeakScores=encyclopediaScorer.getIndividualPeakScores(unitEntry, spectrum, false, uniqueIons.get(peptideModSeq));
					float uniqueScore=score(uniquePeakScores);
					System.out.println(peptideModSeq+": "+uniqueIons.get(peptideModSeq).length+"/"+unitEntry.getIonCount());
					for (int i=0; i<uniquePeakScores.length; i++) {
						if (uniquePeakScores[i]!=null&&uniquePeakScores[i].getScore()>0.0f) {
							System.out.println("\t"+uniquePeakScores[i].getScore()+" for "+uniquePeakScores[i].getTarget()+" ("+uniquePeakScores[i].getTargetMass()+" m/z)");
						}
					}
					String uniqueString=uniqueScore>0?(peptideModSeq+"("+Math.round(uniqueScore)+")"):unknownKey;
					if (bestScore<score) {
						bestScore=score;
						bestTimes.clear();
						bestTimes.put(spectrum.getScanStartTime(), uniqueString);
					} else if (bestScore==score) {
						bestTimes.put(spectrum.getScanStartTime(), uniqueString);
					}
				}
			}

			bestTimes.forEachEntry(new TFloatObjectProcedure<String>() {
				@Override
				public boolean execute(float a, String b) {
					String prev=allBestTimes.get(a);
					if (prev!=null) {
						if (unknownKey.equals(b)) {
							// ignore
						} else if (unknownKey.equals(prev)) {
							allBestTimes.put(a, b);
						} else {
							allBestTimes.put(a, prev+","+b);
						}
					} else {
						allBestTimes.put(a, b);
					}
					return true;
				}
			});

			EValueCalculator calculator=new EValueCalculator(rtScoreMap);
			// System.out.println(peptideModSeq+"\t"+calculator.getMaxRT()+"\t"+calculator.getNegLog10EValue()+"\t"+calculator.getMaxRawScore()); //FIXME
			list.add(calculator.getMaxRT());
		}

		if (allBestTimes.size()>1) {
			float[] times=allBestTimes.keys();
			Arrays.sort(times);
			float range=times[times.length-1]-times[0];
			System.out.println("RANGE: "+range);
			for (int i=0; i<times.length; i++) {
				System.out.println(times[i]/60f+"\t"+allBestTimes.get(times[i]));
			}
			if (range>=60) {
				// System.out.println("multiple\t"+peptideModSeqs.get(0));
				return true;
			}
		} else if (allBestTimes.size()==1) {
			// System.out.println("single\t"+peptideModSeqs.get(0));
		}
		return false;
	}

	private float score(PeakScores[] individualPeakScores) {
		float score=0.0f;
		for (int i=0; i<individualPeakScores.length; i++) {
			if (individualPeakScores[i]!=null&&individualPeakScores[i].getScore()>0.0f) {
				score++;
			}
		}
		return score;
	}
	
	public ArrayList<Stripe> getScanSubset(float minRT, float maxRT) {
		ArrayList<Stripe> subset=new ArrayList<Stripe>();
		for (Stripe scan : stripes) {
			if (scan.getScanStartTime()>=minRT&&scan.getScanStartTime()<=maxRT) {
				subset.add(scan);
			}
		}
		return subset;
	}
}
