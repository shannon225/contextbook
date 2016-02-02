package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;

public class ProphetMixtureModel {
	private Distribution positive;
	private Distribution negative;
	private final boolean fixedMeans;
	
	public ProphetMixtureModel(Distribution positive, Distribution negative, boolean fixedMeans) {
		this.positive=positive;
		this.negative=negative;
		this.fixedMeans=fixedMeans;
	}
	
	public float getProbability(float x) {
		double p=positive.getProbability(x);
		double n=negative.getProbability(x);
		double sum=p+n;
		if (sum==0.0) return 0.0f;
		return (float)(p/sum);
	}
	
	public Distribution getPositive() {
		return positive;
	}
	
	public Distribution getNegative() {
		return negative;
	}
	
	public void train(float[] data, int numIterations) {
		for (int i=0; i<numIterations; i++) {
			runIteration(data);
		}
	}
	
	void runIteration(float[] data) {
		float[] positiveProb=new float[data.length];
		float[] negativeProb=new float[data.length];
		for (int i=0; i<data.length; i++) {
			double positiveValue=positive.getProbability(data[i]);
			double negativeValue=negative.getProbability(data[i]);
			double sum=positiveValue+negativeValue;
			if (sum==0.0) {
				positiveProb[i]=0.0f;
				negativeProb[i]=1.0f;
			} else {
				positiveProb[i]=(float)(positiveValue/sum);
				negativeProb[i]=(float)(negativeValue/sum);
			}
		}
		
		positive=getNewDistribution(positive, data, positiveProb, fixedMeans);
		negative=getNewDistribution(negative, data, negativeProb, fixedMeans);
	}
	
	static Distribution getNewDistribution(Distribution seed, float[] data, float[] weights, boolean fixedMeans) {
		float prior=General.sum(weights);
		//float mean=weightedMean(data, weights);
		//float stdev=weightedStdev(data, weights, mean);
		float mean=trimmedMean(data, weights);
		float stdev=trimmedStdev(data, weights, mean);
		
		if (fixedMeans) {
			return seed.clone(seed.getMean(), stdev, prior);
		} else {
			return seed.clone(mean, stdev, prior);
		}
	}
	
	static float weightedMean(float[] data, float[] weights) {
		return General.sum(data)/General.sum(weights);
	}

	/**
	 * http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
	 * @return
	 */
	static float weightedStdev(float[] data, float[] weights, float mean) {
		float deltaSum=0.0f;
		float weightSum=0.0f;
		int numNonZero=0;
		for (int i=0; i<data.length; i++) {
			if (weights[i]>0.0f) {
				numNonZero++;
				weightSum+=weights[i];
				float delta=data[i]-mean;
				deltaSum+=weights[i]*delta*delta;
			}
		}
		return (float)Math.sqrt((deltaSum)/((numNonZero-1)*weightSum/numNonZero));
	}
	
	static float trimmedMean(float[] data, float[] weights) {
		float sum=0.0f;
		int count=0;
		for (int i=0; i<weights.length; i++) {
			if (weights[i]>=0.5f) {
				sum+=data[i];
				count++;
			}
		}
		return sum/count;
	}

	static float trimmedStdev(float[] data, float[] weights, float mean) {
		float sum=0.0f;
		int count=0;
		for (int i=0; i<weights.length; i++) {
			if (weights[i]>=0.5f) {
				float delta=data[i]-mean;
				sum+=delta*delta;
				count++;
			}
		}
		return (float)Math.sqrt(sum/(count-1));
	}
}
