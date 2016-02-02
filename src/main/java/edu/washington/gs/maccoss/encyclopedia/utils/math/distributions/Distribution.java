package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

public interface Distribution {

	Distribution clone(double mean, double stdev, double prior);

	double getProbability(double x);
	
	double getPDF(double x);

	double getCDF(double x);

	double getMean();

	double getStdev();

	double getPrior();
}