package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;

public class SimpleMixtureModel implements RTProbabilityModel {
	private Distribution positive;
	private Distribution negative;
	
	public SimpleMixtureModel(Distribution positive, Distribution negative) {
		this.positive=positive;
		this.negative=negative;
	}

	@Override
	public float getProbability(float retentionTime, float delta) {
		double pos=positive.getProbability(delta);
		if (pos<0.001f) return 0.0f;
		
		double neg=negative.getProbability(delta);
		double denom = pos+neg;
		if (denom==0) return 0.0f;
		return (float)(pos/denom);
	}
	
	public Distribution getPositive() {
		return positive;
	}
	
	public Distribution getNegative() {
		return negative;
	}
}
