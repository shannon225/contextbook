package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public class GammaDistriution implements Distribution {
	private final double prior;
	

    private static final double halfLog2Pi = 0.5 * Math.log(2.0 * Math.PI);
	private static final double[] lancozsApproximationConstants=new double[] {1.00000000019001, 76.1800917294714,
			-86.5053203294167, 24.0140982408309, -1.23173957245015, 0.00120865097387, -0.000005395239384953};

	private double mean=0.0f;
	private double stdev=0.0f;
	private double minimumValue=0.0f;

	private double normalizedMean=0.0f;
	private double alpha=0.0f;
	private double beta=0.0f;
	private double gammaOfAlpha=0.0f;


	public GammaDistriution(double mean, double stdev, double prior) {
		this.stdev=stdev;
		this.mean=mean;
		this.prior=prior;
		
		normalizedMean=mean-minimumValue;
		alpha=(normalizedMean*normalizedMean)/(stdev*stdev);
		beta=(stdev*stdev)/normalizedMean;
		if (alpha==0.0f) {
			alpha=Double.MIN_VALUE;
		}
		if (beta==0.0f) {
			beta=Double.MIN_VALUE;
		}

		gammaOfAlpha=gamma(alpha);
	}

	public static double gamma(double value) {
		double firstTerm=(Math.sqrt(2*Math.PI))/value;
		double sumOfPs=lancozsApproximationConstants[0];
		for (int i=1; i<lancozsApproximationConstants.length; i++) {
			sumOfPs+=lancozsApproximationConstants[i]/(value+i);
		}

		double lastTerm=(Math.pow(value+5.5, value+0.5)*Math.exp(-(value+5.5)));
		

		return firstTerm*sumOfPs*lastTerm;
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
		return new GammaDistriution(mean, stdev, prior);
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution#getPDF(double)
	 */
	@Override
	public double getPDF(double x) {
		double normalizedValue=x-minimumValue;
		if (normalizedValue<=0) {
			// off the charts in the negative direction
			return 0.0f;
		}

		double powerTerm=0.0f;
		if (normalizedValue!=0.0f) {
			powerTerm=Math.pow(normalizedValue, alpha-1.0f);
		}

		double probability=((powerTerm*Math.exp(-normalizedValue/beta))/(Math.pow(beta, alpha)*gammaOfAlpha));

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
