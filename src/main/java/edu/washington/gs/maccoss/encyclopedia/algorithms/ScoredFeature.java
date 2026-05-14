package edu.washington.gs.maccoss.encyclopedia.algorithms;


public class ScoredFeature {
	private double mz;
	private boolean isDecoy;
	private float primary;
	private float retentionTime;
	private String sequence;
	private String protein;
	private String originalLine;
	
	public ScoredFeature(double mz, boolean isDecoy, float primary, float retentionTime, String sequence, String protein, String originalLine) {
		this.mz=mz;
		this.isDecoy=isDecoy;
		this.primary=primary;
		this.retentionTime=retentionTime;
		this.sequence=sequence;
		this.protein=protein;
		this.originalLine=originalLine;
	}
	
	public double getMz() {
		return mz;
}
	public boolean isDecoy() {
		return isDecoy;
	}
	public float getPrimary() {
		return primary;
	}
	public float getRetentionTime() {
		return retentionTime;
	}
	public String getSequence() {
		return sequence;
	}
	public String getProtein() {
		return protein;
	}
	public String getOriginalLine() {
		return originalLine;
	}
}