package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashMap;
import java.util.Map.Entry;

import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class ScanRangeTracker {
	HashMap<Range, TFloatArrayList> precursorRTsInSecs=new HashMap<Range, TFloatArrayList>();
	HashMap<Range, TFloatArrayList> stripeRTsInSecs=new HashMap<Range, TFloatArrayList>();
	
	public HashMap<Range, Float> getDutyCycleMap() {
		HashMap<Range, Float> dutyCycleMap=new HashMap<Range, Float>();
		for (Entry<Range, TFloatArrayList> entry : stripeRTsInSecs.entrySet()) {
			Range range=entry.getKey();
			TFloatArrayList rts=entry.getValue();
			float[] deltas=General.firstDerivative(rts.toArray());
			float averageDutyCycle=General.mean(deltas);
			dutyCycleMap.put(range, averageDutyCycle);
		}
		return dutyCycleMap;
	}
	
	public HashMap<Range, TFloatArrayList> getStripeRTsInSecs() {
		return stripeRTsInSecs;
	}
	
	public HashMap<Range, TFloatArrayList> getPrecursorRTsInSecs() {
		return precursorRTsInSecs;
	}
	
	/**
	 * NOT THREAD SAFE!
	 * @param range
	 * @param rtInSecs
	 * @return returns if should continue
	 */
	public boolean addRange(Range range, float rtInSecs) {
		TFloatArrayList data=stripeRTsInSecs.get(range);
		if (data==null) {
			data=new TFloatArrayList();
			stripeRTsInSecs.put(range, data);
		}
		
		if (data.size()>1) {
			// already got enough data
			return false;
		}
		data.add(rtInSecs);
		return true;
	}
	
	/**
	 * NOT THREAD SAFE!
	 * @param range
	 * @param rtInSecs
	 */
	public void addPrecursor(Range range, float rtInSecs) {
		TFloatArrayList data=precursorRTsInSecs.get(range);
		if (data==null) {
			data=new TFloatArrayList();
			precursorRTsInSecs.put(range, data);
		}
		
		data.add(rtInSecs);
	}
}
