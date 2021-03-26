package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;

import gnu.trove.list.array.TDoubleArrayList;

public class NearestValueMap<V> {
	TDoubleArrayList keys=new TDoubleArrayList();
	ArrayList<V> values=new ArrayList<>();

	public NearestValueMap() {
	}
	
	public void put(double key, V value) {
		int keyIndex=keys.binarySearch(key);
		// already here
		if (key>=0) values.set(keyIndex, value);
		
		// insertion point
		keyIndex=-(keyIndex+1);
		keys.insert(keyIndex, key);
		values.add(keyIndex, value);
	}
	
	public V get(double key) {
		
		int keyIndex=keys.binarySearch(key);
		if (keyIndex>=0) return values.get(keyIndex);
		

		// insertion point
		keyIndex=-(keyIndex+1);
		
		double deltaUp=Double.MAX_VALUE;
		if (keyIndex>0) {
			deltaUp=Math.abs(keys.get(keyIndex-1)-key);
		}
		double deltaDown=Double.MAX_VALUE;
		if (keyIndex<keys.size()) {
			deltaDown=Math.abs(keys.get(keyIndex)-key);
		}
		
		if (Double.MAX_VALUE==deltaDown&&Double.MAX_VALUE==deltaUp) return null;
		
		if (deltaUp>deltaDown) {
			return values.get(keyIndex-1);
		} else {
			return values.get(keyIndex);
		}
	}
}
