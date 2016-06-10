package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PhosphoLocalizationData {
	private final HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> traces;
	private final HashMap<String, XYTrace[]> uniqueFragmentIons;
	private final HashMap<String, XYPoint> localizationScores;

	public PhosphoLocalizationData(HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> traces, HashMap<String, XYTrace[]> uniqueFragmentIons, HashMap<String, XYPoint> localizationScores) {
		this.traces=traces;
		this.uniqueFragmentIons=uniqueFragmentIons;
		this.localizationScores=localizationScores;
	}
	
	public HashMap<String, XYPoint> getLocalizationScores() {
		return localizationScores;
	}
	
	public HashMap<String, XYTrace[]> getUniqueFragmentIons() {
		return uniqueFragmentIons;
	}
	
	public HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> getTraces() {
		return traces;
	}
	
	public boolean isEmpty() {
		return traces==null||traces.size()==0;
	}
}
