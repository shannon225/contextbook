package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Arrays;

import gnu.trove.impl.Constants;

/** a dynamic, "infinite-length" float array, modeled after TFloatArrayList but with no fixed size **/
public class TFloatDataStore {
    protected float[] _data;
    protected final float no_entry_value;

    public TFloatDataStore() {
    		this(Constants.DEFAULT_CAPACITY, Constants.DEFAULT_FLOAT_NO_ENTRY_VALUE);
    }
    
    public TFloatDataStore(int capacity) {
		this(capacity, Constants.DEFAULT_FLOAT_NO_ENTRY_VALUE);
    }
    
	public TFloatDataStore(float no_entry_value) {
		this(Constants.DEFAULT_CAPACITY, no_entry_value);
	}
    
	public TFloatDataStore(int capacity, float no_entry_value) {
		this._data=new float[capacity];
		this.no_entry_value=no_entry_value;
	}
    
	public void set(int index, float value) {
		ensureCapacity(index+1);
		_data[index]=value;
	}
	
	public float get(int index) {
		if (index>= _data.length) {
			return no_entry_value;
		} else {
			return _data[index];
		}
	}
	
	public float increment(int index, float value) {
		ensureCapacity(index+1);
		if (no_entry_value==_data[index]) {
			_data[index]=value;
		} else {
			_data[index]=_data[index]+value;
		}
		return _data[index];
	}
	
	public float[] toArray(int length) {
		float[] ret=new float[length];
		if (no_entry_value!=0.0f) {
			Arrays.fill(ret, no_entry_value);
		}
		System.arraycopy(_data, 0, ret, 0, Math.min(_data.length, length));
		return ret;
	}
	
	public void replace(float[] data) {
		_data=data.clone();
	}

    protected void ensureCapacity(int capacity) {
        if ( capacity > _data.length ) {
            int newCap = Math.max( _data.length << 1, capacity ); // left shift by 1 doubles capacity
            float[] tmp = new float[ newCap ];
            System.arraycopy( _data, 0, tmp, 0, _data.length );
            _data = tmp;
        }
    }
}
