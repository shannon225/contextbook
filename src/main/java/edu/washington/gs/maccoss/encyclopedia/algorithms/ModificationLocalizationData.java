package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;

public class ModificationLocalizationData {
	private final AmbiguousPeptideModSeq localizationPeptideModSeq;
	private final float localizationScore;
	private final int numberOfMods;
	private final boolean isSiteSpecific;

	public ModificationLocalizationData(AmbiguousPeptideModSeq localizationPeptideModSeq, float localizationScore, int numberOfMods,
			boolean isSiteSpecific) {
		this.localizationPeptideModSeq=localizationPeptideModSeq;
		this.localizationScore=localizationScore;
		this.numberOfMods=numberOfMods;
		this.isSiteSpecific=isSiteSpecific;
	}

	public AmbiguousPeptideModSeq getLocalizationPeptideModSeq() {
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
