package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia.allelespecific;

public class AlleleVariant implements Comparable<AlleleVariant> {
	public static final AlleleVariant EMPTY_VARIANT=new AlleleVariant(0, 0, "", "");
	
	private final int startSite;
	private final int stopSite;
	private final String originalSequence;
	private final String newSequence;

	public AlleleVariant(int site, double modification, char aminoAcid) {
		this.startSite=site;
		this.stopSite=site;
		this.originalSequence=Character.toString(aminoAcid);
		this.newSequence=Character.toString(aminoAcid)+"["+modification+"]";
	}

	public AlleleVariant(int site, char originalSequence, char newSequence) {
		this.startSite=site;
		this.stopSite=site;
		this.originalSequence=Character.toString(originalSequence);
		this.newSequence=Character.toString(newSequence);
	}

	public AlleleVariant(int startSite, int stopSite, String originalSequence, String newSequence) {
		this.startSite=startSite;
		this.stopSite=stopSite;
		this.originalSequence=originalSequence;
		this.newSequence=newSequence;
	}

	public int getStartSite() {
		return startSite;
	}

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
	public int compareTo(AlleleVariant variant ) {
		return this.getStartSite()-variant.getStartSite();
	}
	
}
