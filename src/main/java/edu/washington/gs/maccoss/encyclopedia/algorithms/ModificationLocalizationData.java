package edu.washington.gs.maccoss.encyclopedia.algorithms;

public class ModificationLocalizationData {
	private final String localizationPeptideModSeq;
	private final float localizationScore;
	private final int numberOfMods;
	private final boolean isSiteSpecific;

	public ModificationLocalizationData(String localizationPeptideModSeq, float localizationScore, int numberOfMods,
			boolean isSiteSpecific) {
		this.localizationPeptideModSeq=localizationPeptideModSeq;
		this.localizationScore=localizationScore;
		this.numberOfMods=numberOfMods;
		this.isSiteSpecific=isSiteSpecific;
	}

	public String getLocalizationPeptideModSeq() {
		return localizationPeptideModSeq;
	}

	public float getLocalizationScore() {
		return localizationScore;
	}

	public int getNumberOfMods() {
		return numberOfMods;
	}

	public boolean isSiteSpecific() {
		return isSiteSpecific;
	}

}
