package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

public class UnitDistribution implements Distribution {

	private final double mean;
	private final double stdev;
	private final double prior;
	private final double min;
	private final double max;
	public UnitDistribution(double mean, double stdev, double prior, double min, double max) {
		this.mean=mean;
		this.stdev=stdev;
		this.prior=prior;
		this.min=min;
		this.max=max;
	}
	
	@Override
	public String getName() {
		return "Unit";
	}
	
	@Override
	public String toString() {
		return "Unit, m:"+mean+", sd:"+stdev+", p:"+prior;
	}
	
	public double getMean() {
		return mean;
	}
	public double getStdev() {
		return stdev;
	}
	public double getPrior() {
		return prior;
	}
	@Override
	public Distribution clone(double mean, double stdev, double prior) {
		return new UnitDistribution(mean, stdev, prior, min, max);
	}
	@Override
	public double getProbability(double x) {
		return getPDF(x)*getPrior();
	}
	@Override
	public double getPDF(double x) {
		return 1.0/(max-min);
	}
	@Override
	public double getCDF(double x) {
		return (x-min)/(max-min);
	}
	
	
}
