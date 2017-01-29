package edu.washington.gs.maccoss.encyclopedia.utils;

import gnu.trove.map.hash.TIntFloatHashMap;

public class SparseIndexMap extends TIntFloatHashMap {

    public SparseIndexMap() {
		super();
	}

	public SparseIndexMap(int initialCapacity) {
		super(initialCapacity);
	}
	
	public void multiplyAllValues(float value) {
		for (int i=0; i<_values.length; i++) {
			_values[i]=_values[i]*value;
		}
	}

	public float putIfGreater( int key, float value ) {
        int index = insertKey( key );
        if (index < 0) {
        	if (value<=_values[-index - 1]) {
                return _values[-index - 1];
        	}
        }
        return doPut( key, value, index );
    }
    
    // copy of private method
    private float doPut( int key, float value, int index ) {
        float previous = no_entry_value;
        boolean isNewMapping = true;
        if ( index < 0 ) {
            index = -index -1;
            previous = _values[index];
            isNewMapping = false;
        }
        _values[index] = value;

        if (isNewMapping) {
            postInsertHook( consumeFreeSlot );
        }

        return previous;
    }
}
