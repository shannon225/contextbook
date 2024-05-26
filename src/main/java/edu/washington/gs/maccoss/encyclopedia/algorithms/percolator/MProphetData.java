package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

public class MProphetData implements Comparable<MProphetData> {
	private final String id;
	private final String sequence;
	private final String protein;
	private final float[] data;
	private final boolean isDecoy;
	
	public MProphetData(String id, String sequence, String protein, float[] data, boolean isDecoy) {
		this.id = id;
		this.sequence = sequence;
		this.protein = protein;
		this.data = data;
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
