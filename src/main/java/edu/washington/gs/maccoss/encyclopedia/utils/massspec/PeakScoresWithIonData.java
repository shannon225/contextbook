package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Optional;

public class PeakScoresWithIonData extends PeakScores {
	private final Optional<Float> ionMobility;

	public PeakScoresWithIonData(float score, FragmentIon target, float deltaMass, Optional<Float> ionMobility) {
		super(score, target, deltaMass);
		this.ionMobility = ionMobility;
	}
	
	public Optional<Float> getIonMobility() {
		return ionMobility;
	}
}
