package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

public class CosineGaussian implements Distribution {
	private final double prior;
	private final double stdev;
	private final double mean;

	private final double A;
	private final double M;
	private final double min;
	private final double max;
	
	public CosineGaussian(double mean, double stdev, double prior) {
		this.stdev=stdev;
		this.mean=mean;
		this.prior=prior;

		A=Math.PI/(8.0*stdev);
		M=Math.PI/(4.0*stdev);
		this.min=mean-2.0*stdev;
		this.max=mean+2.0*stdev;
	}
	
	@Override
	public String getName() {
		return "Cosine Distribution";
	}
	
	@Override
	public String toString() {
		return "Cosine, m:"+mean+", sd:"+stdev+", p:"+prior;
	}
	
	@Override
	public Distribution clone(double mean, double stdev, double prior) {
		return new CosineGaussian(mean, stdev, prior);
	}

	@Override
	public double getProbability(double x) {
		return getPDF(x)*getPrior();
	}

	@Override
	public double getPDF(double x) {
		if (x<min) return 0.0f;
		if (x>max) return 0.0f;		
		return A*cosine((mean-x)*M);
	}

	/**
	 * CDF is estimated using logistic approximation:
	 * Bowling et al 2009, doi:10.3926/jiem.2009.v2n1.p114-127
	 */
	@Override
	public double getCDF(double x) {
		if (x<min) return 0.0f;
		if (x>max) return 1.0f;
		
		double d=(x-mean)/stdev;
		return 1.0/(1.0+Math.exp(-0.07056*d*d*d-1.5976*d));
	}

	@Override
	public double getMean() {
		return mean;
	}

	@Override
	public double getStdev() {
		return stdev;
	}

	@Override
	public double getPrior() {
		return prior;
	}

	/**
	 * fast cosine detailed here:
	 * http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648
	 */
	private static final double pi=Math.PI;
	private static final double halfPi=pi/2.0;
	private static final double B=4.0/pi;
	private static final double C=-4.0/(pi*pi);
	private static final double P=0.225;

	public static double cosine(double v) {
		double x=v+halfPi;
		double y=B*x+C*x*Math.abs(x);
		y=P*(y*Math.abs(y)-y)+y;
		return y;
	}
}
