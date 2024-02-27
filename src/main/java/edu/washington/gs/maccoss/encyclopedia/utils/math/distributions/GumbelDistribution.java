package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class GumbelDistribution implements Distribution {
	private final double prior;
	private static final float eulersConstant=0.57721566490153286f;

	private double mean=0.0f;
	private double stdev=0.0f;

	private double alpha=0.0f;
	private double beta=0.0f;

	public GumbelDistribution(double mean, double stdev, double prior) {
		this.stdev=stdev;
		this.mean=mean;
		this.prior=prior;

		this.beta=Math.sqrt(6.0f)/Math.PI*stdev;
		this.alpha=mean-beta*eulersConstant;

		if (beta==0.0f) {
			beta=Double.MIN_VALUE;
		}

	}
	
	@Override
	public String getName() {
		return "Gaussian";
	}
	
	@Override
	public String toString() {
		return "Gaussian, m:"+mean+", sd:"+stdev+", p:"+prior;
	}

	@Override
	public double getPrior() {
		return prior;
	}
	
	@Override
	public double getProbability(double x) {
		return getPDF(x)*getPrior();
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#clone(double, double)
	 */
	@Override
	public Distribution clone(double mean, double stdev, double prior) {
		return new GumbelDistribution(mean, stdev, prior);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getPDF(double)
	 */
	@Override
	public double getPDF(double x) {
		double internal=(alpha-x)/beta;
		double probability=Math.exp(internal-Math.exp(internal))/beta;

		if (Double.isNaN(probability)) {
			return 0.0f;
		} else if (Double.isInfinite(probability)){
			return 1.0f;
		} else {
			return (float)probability;
		}
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getCDF(double)
	 */
	@Override
	public double getCDF(double x) {
		throw new EncyclopediaException("CDF not implemented!");
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getMean()
	 */
	@Override
	public double getMean() {
		return mean;
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getStdev()
	 */
	@Override
	public double getStdev() {
		return stdev;
	}
}
