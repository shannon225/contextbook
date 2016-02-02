package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

public class Gaussian implements Distribution {
	private final double prior;
	private final double stdev;
	private final double mean;
	private final double constant;
	private final double doubleVariance;

	public Gaussian(double mean, double stdev, double prior) {
		this.stdev=stdev;
		this.mean=mean;
		this.prior=prior;
		
		double variance=stdev*stdev;
		this.constant=1.0/(stdev*Math.sqrt(2*Math.PI));
		this.doubleVariance=2.0*variance;
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
		return new Gaussian(mean, stdev, prior);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getPDF(double)
	 */
	@Override
	public double getPDF(double x) {
		double delta=x-mean;
		return constant*Math.exp(-(delta*delta)/doubleVariance);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getCDF(double)
	 */
	@Override
	public double getCDF(double x) {
		double sign=1;
		if (x<0) sign=-1;

		double result=0.5*(1.0+sign*erf(Math.abs(x)/sqrt2));
		return result;
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

	// A&S formula 7.1.26
	private static final double a1=0.254829592;
	private static double a2=-0.284496736;
	private static double a3=1.421413741;
	private static double a4=-1.453152027;
	private static double a5=1.061405429;
	private static double p=0.3275911;
	private static double sqrt2=Math.sqrt(2.0);

	private static double erf(double x) {
		x=Math.abs(x);
		double t=1/(1+p*x);
		
		// Horner's method, takes O(n) operations for nth order polynomial
		return 1-((((((a5*t+a4)*t)+a3)*t+a2)*t)+a1)*t*Math.exp(-1*x*x);
	}
}
