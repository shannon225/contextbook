package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class Gaussian {
	private final double stdev;
	private final double mean;
	private final double constant;
	private final double doubleVariance;

	public Gaussian(double mean, double stdev) {
		this.stdev=stdev;
		this.mean=mean;
		double variance=stdev*stdev;
		this.constant=1.0/(stdev*Math.sqrt(2*Math.PI));
		this.doubleVariance=2.0*variance;
	}

	public double getY(double x) {
		double delta=x-mean;
		return constant*Math.exp(-(delta*delta)/doubleVariance);
	}
	
	public double getMean() {
		return mean;
	}
	
	public double getStdev() {
		return stdev;
	}
}
