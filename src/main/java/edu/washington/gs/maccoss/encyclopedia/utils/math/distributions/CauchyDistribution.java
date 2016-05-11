package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

public class CauchyDistribution implements Distribution {
	private final double mean;
	private final double gamma;
	private final double prior;
	

	public CauchyDistribution(double mean, double gamma, double prior) {
		this.mean=mean;
		this.gamma=gamma;
		this.prior=prior;
	}

	@Override
	public String getName() {
		return "Cauchy";
	}

	@Override
	public Distribution clone(double mean, double stdev, double prior) {
		return new CauchyDistribution(mean, gamma, prior);
	}

	@Override
	public double getProbability(double x) {
		return getPDF(x)*getPrior();
	}

	@Override
	public double getPDF(double x) {
		double inner=(x-mean)/gamma;
		return 1.0/(Math.PI*gamma*(1+inner*inner));
	}

	@Override
	public double getCDF(double x) {
		return 1.0/Math.PI*Math.atan(x/gamma)+0.5;
	}

	@Override
	public double getMean() {
		return mean;
	}

	@Override
	public double getStdev() {
		return gamma*2.0;
	}

	@Override
	public double getPrior() {
		return prior;
	}

}
