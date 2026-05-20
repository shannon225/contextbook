package edu.washington.gs.maccoss.encyclopedia.context;


public class ScoredFeature {
	private double mz;
	private byte charge;
	private boolean isDecoy;
	private float primary;
	private float retentionTime;
	private String sequence;
	private String protein;
	private String originalLine;
	private boolean isBackground;

	public ScoredFeature(double mz, boolean isDecoy, float primary, float retentionTime, String sequence, String protein, String originalLine) {
		this.mz=mz;
		this.isDecoy=isDecoy;
		this.primary=primary;
		this.retentionTime=retentionTime;
		this.sequence=sequence;
		this.protein=protein;
		this.originalLine=originalLine;
	}

	public ScoredFeature(double mz, boolean isDecoy, float primary, float retentionTime, String sequence, String protein, String originalLine, boolean isBackground) {
		this.mz=mz;
		this.isDecoy=isDecoy;
		this.primary=primary;
		this.retentionTime=retentionTime;
		this.sequence=sequence;
		this.protein=protein;
		this.originalLine=originalLine;
		this.isBackground=isBackground;

	}

<<<<<<< HEAD
<<<<<<< HEAD
	public ScoredFeature(double mz, byte charge, boolean isDecoy, float primary, float retentionTime, String sequence, String protein, String originalLine) {
=======
	public ScoredFeature(double mz, byte charge, boolean isDecoy, float primary, float retentionTime, String sequence, String protein, String originalLine, boolean isBackground) {
>>>>>>> f44678a1 (Added a class that will process features with Encyclopedia without running Percolator, and export them as pin.tsv files.)
=======
	public ScoredFeature(double mz, byte charge, boolean isDecoy, float primary, float retentionTime, String sequence, String protein, String originalLine) {
>>>>>>> 8b8d896c (Added sequences to the mass lists for the TargetedBoostrapper class.)
		this.mz=mz;
		this.charge=charge;
		this.isDecoy=isDecoy;
		this.primary=primary;
		this.retentionTime=retentionTime;
		this.sequence=sequence;
		this.protein=protein;
		this.originalLine=originalLine;
<<<<<<< HEAD
<<<<<<< HEAD
//		this.isBackground=isBackground; 
=======
		this.isBackground=isBackground; 
>>>>>>> f44678a1 (Added a class that will process features with Encyclopedia without running Percolator, and export them as pin.tsv files.)
=======
//		this.isBackground=isBackground; 
>>>>>>> 8b8d896c (Added sequences to the mass lists for the TargetedBoostrapper class.)
	}


	public double getMz() {
		return mz;
	}

	public byte getCharge() {
		return charge;
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
	public boolean isBackground() {
		return isBackground;
	}
}