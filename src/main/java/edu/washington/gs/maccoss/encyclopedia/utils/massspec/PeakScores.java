package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public class PeakScores {
	private final float score;
	private final double targetMass;
	private final float deltaMass;
	public PeakScores(float score, double targetMass, float deltaMass) {
		this.score = score;
		this.targetMass=targetMass;
		this.deltaMass = deltaMass;
	}
	public double getTargetMass() {
		return targetMass;
	}
	
	public float getDeltaMass() {
		return deltaMass;
	}
	
	public float getScore() {
		return score;
	}
	
	public static float sumScores(PeakScores[] scores) {
		float score=0.0f;
		for (int i = 0; i < scores.length; i++) {
			if (scores[i]!=null) {
				score+=scores[i].score;
			}
		}
		return score;
	}
}
