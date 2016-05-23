package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.special.Beta;

import edu.washington.gs.maccoss.encyclopedia.algorithms.DotProduct;
import edu.washington.gs.maccoss.encyclopedia.algorithms.EValueCalculator;
import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PeptideQuantExtractorTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaOneScorer;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.math.BinomialCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.map.hash.TFloatObjectHashMap;
import gnu.trove.procedure.TFloatObjectProcedure;

public class PhosphoLocalizer {
	private static final String unknownKey="unknown";
	private final String filename;
	private final SearchParameters params;
	private final PSMData psmdata;
	private final ArrayList<Stripe> stripes;
	private final EncyclopediaOneScorer encyclopediaScorer;
	private final DotProduct dotProduct;

	public PhosphoLocalizer(String filename, SearchParameters params, PSMData psmdata, ArrayList<Stripe> stripes) {
		this.filename=filename;
		this.params=params;
		this.psmdata=psmdata;
		this.stripes=stripes;
		encyclopediaScorer=new EncyclopediaOneScorer(params, null); // not using aux scoring
		dotProduct=new DotProduct(params.getFragmentTolerance());
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
		MassTolerance fragmentTolerance=params.getFragmentTolerance();
		float prior=(float)fragmentTolerance.getTolerance(500.0)/100.0f; // equivalent to high mass accuracy a-score
		
		float duration=6*60f; // search for 6 minutes

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
		
		for (Entry<String, FragmentationModel> entry : entryMap.entrySet()) {
			String peptideModSeq=entry.getKey();
			LibraryEntry unitEntry=entry.getValue().getUnitSpectrum(filename, new HashSet<String>(), precursorCharge, retentionTime, params);

			TFloatFloatHashMap rtScoreMap=new TFloatFloatHashMap();
			TFloatArrayList scores=new TFloatArrayList();
			for (Stripe spectrum : stripes) {
				float score=encyclopediaScorer.score(unitEntry, spectrum);
				scores.add(score);
				rtScoreMap.put(spectrum.getScanStartTime(), score);
			}
			EValueCalculator calculator=new EValueCalculator(rtScoreMap);

			TFloatFloatHashMap uniqueRtScoreMap=new TFloatFloatHashMap();
			FragmentIon[] targets=uniqueIons.get(peptideModSeq);
			int bestFound=0;
			int bestLookedFor=0;
			float bestEvalue=0.0f;
			for (int i=0; i<stripes.size(); i++) {
				Stripe spectrum=stripes.get(i);
				float score=scores.get(i);
				float evalue=calculator.getNegLog10EValue(score);
				if (evalue>-1.0f) {
					double[] masses=FragmentIon.getMasses(targets);
					float[] unitIntensities=new float[masses.length];
					Arrays.fill(unitIntensities, 1.0f);
					LibraryEntry localUnit=new LibraryEntry(unitEntry.getSource(), unitEntry.getAccessions(), unitEntry.getSpectrumIndex(), precursorMZ, precursorCharge, peptideModSeq, 1, spectrum.getScanStartTime(), score, masses, unitIntensities);
					
					TransitionRefinementData data=PeptideQuantExtractorTask.quantifyPeptide(dotProduct, localUnit, false, stripes);
					if (data==null) continue;
					
					int found=0;
					int lookedFor=0;
					ArrayList<float[]> chromatograms=data.getChromatograms();
					for (float[] fs : chromatograms) {
						for (int j=0; j<fs.length; j++) {
							lookedFor++;
							if (fs[j]>0.0f) {
								found++;
							}
						}
					}
					if (found>bestFound) {
						bestFound=found;
						bestLookedFor=lookedFor;
						bestEvalue=evalue;
					}

					if (lookedFor>0&&found>0) {
						double pValue=cumulativeBinomalGreaterOrEqual(lookedFor, found, prior);
						uniqueRtScoreMap.put(spectrum.getScanStartTime(), -10.0f*(float)Log.log10(pValue));
						System.out.println(peptideModSeq+" --> "+spectrum.getScanStartTime()/60.0f+", "+pValue+" ("+bestFound+"/"+bestLookedFor+")");
					}
				}
			}
			EValueCalculator uniqueCalculator=new EValueCalculator(uniqueRtScoreMap);

			float bestRT=uniqueCalculator.getMaxRT();
			System.out.println("FINAL: "+peptideModSeq+" --> "+bestRT/60.0f+", "+uniqueCalculator.getMaxRawScore()+" ("+bestFound+"/"+bestLookedFor+") evalue: "+bestEvalue);
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
	private static double cumulativeBinomalGreaterOrEqual(int n, int k, double p) {
		BinomialCalculator dist=new BinomialCalculator(n, p);
		return dist.cumulativeProbabilityGreaterThan(k-1);
	}

	public static int countMatchingTargets(Stripe spectrum, FragmentIon[] targets, MassTolerance tolerance) {
		double[] acquiredMasses=spectrum.getMassArray();
		float[] acquiredIntensities=spectrum.getIntensityArray();
		
		int count=0;
		for (FragmentIon fragmentIon : targets) {
			double target=fragmentIon.mass;
			int[] indicies=tolerance.getIndicies(acquiredMasses, target);
			for (int i=0; i<indicies.length; i++) {
				if (acquiredIntensities[indicies[i]]>0.0f) {
					count++;
				}
			}
		}
		return count;
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
