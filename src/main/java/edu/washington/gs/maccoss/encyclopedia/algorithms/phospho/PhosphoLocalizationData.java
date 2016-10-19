package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PhosphoLocalizationData {
	private final HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> scoreTraces;
	private final HashMap<String, HashMap<FragmentIon, XYTrace>> uniqueFragmentIons;
	private final HashMap<String, HashMap<FragmentIon, XYTrace>> otherFragmentIons;
	private final HashMap<String, FragmentIon[]> uniqueTargetFragments;
	private final HashMap<String, XYPoint> localizationScores;
	private final HashMap<String, TransitionRefinementData> passingForms;

	public PhosphoLocalizationData(HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> scoreTraces, HashMap<String, HashMap<FragmentIon, XYTrace>> uniqueFragmentIons, HashMap<String, HashMap<FragmentIon, XYTrace>> otherFragmentIons, HashMap<String, FragmentIon[]> uniqueTargetFragments, HashMap<String, XYPoint> localizationScores, HashMap<String, TransitionRefinementData> passingForms) {
		this.scoreTraces=scoreTraces;
		this.uniqueFragmentIons=uniqueFragmentIons;
		this.uniqueTargetFragments=uniqueTargetFragments;
		this.otherFragmentIons=otherFragmentIons;
		this.localizationScores=localizationScores;
		this.passingForms=passingForms;
	}
	
	public HashMap<String, XYPoint> getLocalizationScores() {
		return localizationScores;
	}
	
	public HashMap<String, HashMap<FragmentIon, XYTrace>> getUniqueFragmentIons() {
		return uniqueFragmentIons;
	}
	
	public HashMap<String, HashMap<FragmentIon, XYTrace>> getOtherFragmentIons() {
		return otherFragmentIons;
	}
	
	public HashMap<String, FragmentIon[]> getUniqueTargetFragments() {
		return uniqueTargetFragments;
	}
	
	public HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> getScoreTraces() {
		return scoreTraces;
	}
	
	public HashMap<String, TransitionRefinementData> getPassingForms() {
		return passingForms;
	}
	
	public boolean isEmpty() {
		return scoreTraces==null||scoreTraces.size()==0;
	}
}
