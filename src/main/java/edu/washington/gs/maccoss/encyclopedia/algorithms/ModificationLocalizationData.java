package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;

public class ModificationLocalizationData {
	public static final ModificationLocalizationData POISON_RESULT=new ModificationLocalizationData(null, 0.0f, 0.0f, 0, false, false, null, 0.0f, 0.0f);
	
	private final AmbiguousPeptideModSeq localizationPeptideModSeq;
	private final float retentionTimeApexInSeconds;
	private final float localizationScore;
	private final int numberOfMods;
	private final boolean isSiteSpecific;
	private final boolean isLocalized;
	private final FragmentIon[] localizingIons;
	private final float localizingIntensity;
	private final float totalIntensity;

	public ModificationLocalizationData(AmbiguousPeptideModSeq localizationPeptideModSeq, float retentionTimeApexInSeconds, float localizationScore, int numberOfMods, boolean isSiteSpecific, boolean isLocalized,
			FragmentIon[] localizingIons, float localizingIntensity, float totalIntensity) {
		this.localizationPeptideModSeq=localizationPeptideModSeq;
		this.retentionTimeApexInSeconds=retentionTimeApexInSeconds;
		this.localizationScore=localizationScore;
		this.numberOfMods=numberOfMods;
		this.isSiteSpecific=isSiteSpecific;
		this.isLocalized=isLocalized;
		this.localizingIons=localizingIons;
		this.localizingIntensity=localizingIntensity;
		this.totalIntensity=totalIntensity;
	}

	public float getRetentionTimeApexInSeconds() {
		return retentionTimeApexInSeconds;
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
	
	public boolean isLocalized() {
		return isLocalized;
	}

	public FragmentIon[] getLocalizingIons() {
		return localizingIons;
	}

	public float getLocalizingIntensity() {
		return localizingIntensity;
	}

	public float getTotalIntensity() {
		return totalIntensity;
	}
}
