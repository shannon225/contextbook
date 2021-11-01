package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

public class LaplaceDistribution implements Distribution {
	private final double prior;
	private final double stdev;
	private final double mean;
	private final double b;

	public LaplaceDistribution(double mean, double stdev, double prior) {
		this.stdev=stdev;
		this.mean=mean;
		this.prior=prior;
		b=stdev/Math.sqrt(2.0);
	}
	
	@Override
	public String getName() {
		return "Laplace";
	}
	
	@Override
	public String toString() {
		return "Laplace, m:"+mean+", sd:"+stdev+", p:"+prior;
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
		return new LaplaceDistribution(mean, stdev, prior);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getPDF(double)
	 */
	@Override
	public double getPDF(double x) {
		return (0.5/b)*Math.exp(-Math.abs(x-mean)/b);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getCDF(double)
	 */
	@Override
	public double getCDF(double x) {
		if (x>mean) {
			return 1-0.5*Math.exp(-Math.abs(x-mean)/b);
		} else {
			return 0.5*Math.exp(-Math.abs(x-mean)/b);
		}
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
