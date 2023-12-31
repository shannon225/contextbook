package edu.washington.gs.maccoss.encyclopedia.utils;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import gnu.trove.map.hash.TIntObjectHashMap;

@SuppressWarnings("rawtypes")
public class SparseIndexMap extends TIntObjectHashMap {

	public SparseIndexMap() {
		super();
	}

	public SparseIndexMap(int initialCapacity) {
		super(initialCapacity);
	}
	
	public Peak getPeak(int index) {
		return (Peak)get(index);
	}

	public void multiplyAllValues(float value) {
		for (int i=0; i<_values.length; i++) {
			if (_values[i]!=null) {
				Peak peak=(Peak)_values[i];
				_values[i]=new Peak(peak.mass, peak.intensity*value, peak.ionMobility);
			}
		}
	}

	public Peak adjustOrPutValue(int key, double mass, float intensity) {
		return adjustOrPutValue(key, mass, intensity, null);
	}
	public Peak adjustOrPutValue(int key, double mass, float intensity, Float ionmobility) {
		int index=insertKey(key);
		final boolean isNewMapping;
		final Peak oldValue;
		if (index<0) {
			index=-index-1;
			oldValue=(Peak)_values[index];
			_values[index]=new Peak(oldValue.mass, oldValue.intensity+intensity, oldValue.ionMobility);
			isNewMapping=false;
		} else {
			oldValue=(Peak)_values[index];
			_values[index]=new Peak(mass, intensity, ionmobility);
			isNewMapping=true;
		}

		if (isNewMapping) {
			postInsertHook(consumeFreeSlot);
		}

		return oldValue;
	}

	public Peak putIfGreater(int key, double mass, float intensity, Float ionmobility) {
		int index=insertKey(key);
		if (index<0) {
			Peak peak=(Peak)_values[-index-1];
			if (intensity<=peak.intensity) {
				return peak;
			}
		}
		return doPut(key, new Peak(mass, intensity, ionmobility), index);
	}

	// copy of private method
	private Peak doPut(int key, Peak value, int index) {
		Peak previous=null;
		boolean isNewMapping=true;
		if (index<0) {
			index=-index-1;
			previous=(Peak)_values[index];
			isNewMapping=false;
		}
		_values[index]=value;

		if (isNewMapping) {
			postInsertHook(consumeFreeSlot);
		}

		return previous;
	}
}
