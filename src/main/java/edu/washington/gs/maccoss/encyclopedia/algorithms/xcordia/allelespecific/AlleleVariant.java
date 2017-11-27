package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

public class AlleleVariant implements Comparable<AlleleVariant> {
	public static final AlleleVariant EMPTY_VARIANT=new AlleleVariant(0, 0, "", "");
	
	private final int startSite; // count starting from 1!
	private final int stopSite;
	private final String originalSequence;
	private final String newSequence;

	/**
	 * Note: count starting from 1!
	 */
	public AlleleVariant(int site, double modification, char aminoAcid) {
		this.startSite=site;
		this.stopSite=site;
		this.originalSequence=Character.toString(aminoAcid);
		this.newSequence=Character.toString(aminoAcid)+"["+modification+"]";
	}

	/**
	 * Note: count starting from 1!
	 */
	public AlleleVariant(int site, char originalSequence, char newSequence) {
		this.startSite=site;
		this.stopSite=site;
		this.originalSequence=Character.toString(originalSequence);
		this.newSequence=Character.toString(newSequence);
	}

	/**
	 * Note: count starting from 1!
	 */
	public AlleleVariant(int startSite, int stopSite, String originalSequence, String newSequence) {
		this.startSite=startSite;
		this.stopSite=stopSite;
		this.originalSequence=originalSequence;
		this.newSequence=newSequence;
	}
	
	public String toString() {
		return startSite+"("+originalSequence+"/"+newSequence+")";
	}

	/**
	 * Note: count starting from 1!
	 * @return
	 */
	public int getStartSite() {
		return startSite;
	}

	/**
	 * Note: count starting from 1!
	 * @return
	 */
	public int getStopSite() {
		return stopSite;
	}

	public String getOriginalSequence() {
		return originalSequence;
	}

	public String getNewSequence() {
		return newSequence;
	}
	
	@Override
	public int hashCode() {
		return getStartSite()+16807*getStopSite();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==null||!(obj instanceof AlleleVariant)) return false;
		return compareTo((AlleleVariant)obj)==0;
	}

	@Override
	public int compareTo(AlleleVariant variant) {
		if (variant==null) return 1;
		int c=Integer.compare(this.getStartSite(), variant.getStartSite());
		if (c!=0) return c;
		c=Integer.compare(this.getStopSite(), variant.getStopSite());
		if (c!=0) return c;
		c=this.getOriginalSequence().compareTo(variant.getOriginalSequence());
		if (c!=0) return c;
		return this.getNewSequence().compareTo(variant.getNewSequence());
	}
	
}
