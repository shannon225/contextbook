package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public class PeakScores {
	private final float score;
	private final FragmentIon target;
	private final float deltaMass;
	public PeakScores(float score, FragmentIon target, float deltaMass) {
		this.score = score;
		this.target=target;
		this.deltaMass = deltaMass;
	}
	public String toString() {
		return target.toString()+"="+score;
	}
	public double getTargetMass() {
		return target.getMass();
	}
	
	public FragmentIon getTarget() {
		return target;
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
