package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.map.hash.TFloatFloatHashMap;

public class PhosphoLocalizationData {
	private final HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> traces;

	public PhosphoLocalizationData(HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> traces) {
		this.traces=traces;
	}
	
	public HashMap<String, Pair<TFloatFloatHashMap, TFloatFloatHashMap>> getTraces() {
		return traces;
	}
	
	public boolean isEmpty() {
		return traces==null||traces.size()==0;
	}
}
