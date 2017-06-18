package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;

public class LocalizationData {

	private final float bestLocalizationScore;
	private final Stripe apex;
	private final boolean wasLocalized;

	public LocalizationData(float bestLocalizationScore, Stripe apex, boolean wasLocalized) {
		this.bestLocalizationScore=bestLocalizationScore;
		this.apex=apex;
		this.wasLocalized=wasLocalized;
	}

	public Stripe getApex() {
		return apex;
	}
	public float getBestLocalizationScore() {
		return bestLocalizationScore;
	}
	public boolean wasLocalized() {
		return wasLocalized;
	}
}
