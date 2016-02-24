package edu.washington.gs.maccoss.encyclopedia.utils.math.randomforest;

import gnu.trove.list.array.TFloatArrayList;

public class Record implements Comparable<Record> {
	private final String sampleName;
	private final boolean state;

	private final float[] data;

	public Record(String sampleName, boolean state, float[] data) {
		this.sampleName=sampleName;
		this.state=state;
		this.data=data;
	}
	
	@Override
	/**
	 * not really comparable! Just compares hashcode as a cheat!
	 */
	public int compareTo(Record o) {
		if (o==null) return 1;
		return hashCode()-o.hashCode();
	}

	public String getSampleName() {
		return sampleName;
	}

	public float[] getData() {
		return data;
	}

	public float[] getTrimmedData(boolean[] useData) {
		TFloatArrayList list=new TFloatArrayList();
		for (int i=0; i<useData.length; i++) {
			if (useData[i])
				list.add(data[i]);
		}

		return list.toArray();
	}

	public int getDataLength() {
		return data.length;
	}

	public float getDataValue(int index) {
		return data[index];
	}

	public boolean state() {
		return state;
	}
}