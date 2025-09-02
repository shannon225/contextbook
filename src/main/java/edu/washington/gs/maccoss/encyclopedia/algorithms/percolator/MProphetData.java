package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

public class MProphetData implements Comparable<MProphetData> {
	private final String id;
	private final String sequence;
	private final String protein;
	private final float[] data;
	private final boolean isDecoy;
	
	// forces N/A data to be 0.0
	public MProphetData(String id, String sequence, String protein, float[] startingData, boolean isDecoy) {
		this.id = id;
		this.sequence = sequence;
		this.protein = protein;
		this.data = new float[startingData.length];
		for (int i = 0; i < data.length; i++) {
			if (Float.isNaN(startingData[i])) {
				data[i]=0.0f;
			} else {
				data[i]=startingData[i];
			}
		}
		
		this.isDecoy=isDecoy;
	}
	
	@Override
	public int compareTo(MProphetData o) {
		return id.compareTo(o.id);
	}
	@Override
	public int hashCode() {
		return id.hashCode();
	}
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof MProphetData) return compareTo((MProphetData)obj)==0;
		return false;
	}
	
	public float[] getData() {
		return data;
	}
	public String getId() {
		return id;
	}
	public String getProtein() {
		return protein;
	}
	public String getSequence() {
		return sequence;
	}
	public boolean isDecoy() {
		return isDecoy;
	}
}
