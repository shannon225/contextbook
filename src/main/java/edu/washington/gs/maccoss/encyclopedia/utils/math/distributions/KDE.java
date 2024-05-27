package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.WeightedValue;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class KDE implements Distribution {
	private final ArrayList<Distribution> data;
	private final Range range;
	private final int numberOfBins;
	private final double[] histogram;
	private final double sumPriors;
	private final double sumHistogram;
	private final double prior;
	private final double mean;
	private final double stdev;
	private final PolynomialSplineFunction spline;
	
	public KDE(double[] values, double prior) {
		this(values, prior, 100);
	}
	
	public KDE(float[] values, double prior) {
		this(values, prior, 100);
	}
	
	public KDE(double[] values, double prior, int numberOfBins) {
		this(getUnitWeightedValues(values), prior, numberOfBins);
	}
	
	public KDE(float[] values, double prior, int numberOfBins) {
		this(getUnitWeightedValues(General.toDoubleArray(values)), prior, numberOfBins);
	}
	
	private static ArrayList<WeightedValue> getUnitWeightedValues(double[] values) {
		ArrayList<WeightedValue> list=new ArrayList<WeightedValue>();
		for (int i = 0; i < values.length; i++) {
			list.add(new WeightedValue(values[i], 1));
		}
		return list;
	}
	public KDE(ArrayList<WeightedValue> values, double prior) {
		this(values, prior, 100);
	}

	public KDE(ArrayList<WeightedValue> values, double prior, int numberOfBins) {
		this.prior=prior;
		this.numberOfBins=numberOfBins;
		stdev=WeightedValue.stdev(values);
		mean=WeightedValue.mean(values);
		// Silverman's (1986) rule of thumb (wikipedia)
		double bandwidth=stdev*Math.pow(4.0/3.0/values.size(), 1.0/5.0);

		data=new ArrayList<Distribution>();
		for (WeightedValue value : values) {
			// division by 2.3548 converts bandwidth (fwhm) to stdev for
			// gaussians
			data.add(new CosineGaussian(value.getValue(), bandwidth/2.3548, value.getWeight()));
		}

		float min=Float.MAX_VALUE;
		float max=-Float.MAX_VALUE;
		for (Distribution dist : data) {
			float localMin=(float)(dist.getMean()-2.0f*dist.getStdev());
			float localMax=(float)(dist.getMean()+2.0f*dist.getStdev());
			if (min>localMin) min=localMin;
			if (max<localMax) max=localMax;
		}
		range=new Range(min, max);
		double binsize=Math.max(1e-4, (max-min)/numberOfBins); // force increasing

		histogram=new double[numberOfBins];
		double[] binValues=new double[histogram.length];
		for (int i=0; i<histogram.length; i++) {
			binValues[i]=i*binsize+min;
		}

		double total=0.0;
		for (Distribution dist : data) {
			total+=dist.getPrior();
			double localMin=(double)(dist.getMean()-2.0f*dist.getStdev());
			double localMax=(double)(dist.getMean()+2.0f*dist.getStdev());
			int startIndex=Math.max(0, (int)Math.floor((localMin-min)/binsize));
			int stopIndex=Math.min(histogram.length-1, (int)Math.ceil((localMax-min)/binsize));

			for (int i=startIndex; i<=stopIndex; i++) {
				double probability=(double)dist.getProbability(binValues[i]);
				if (!Double.isNaN(probability)&&!Double.isInfinite(probability)) {
					histogram[i]+=probability;
				}
			}
		}
		sumHistogram=General.sum(histogram);
		sumPriors=total;

		if (data.size()>10) {
			spline=new SplineInterpolator().interpolate(binValues, histogram);
		} else {
			spline=null;
		}
	}
	
	@Override
	public String getName() {
		return "KDE";
	}

	public int getBin(double value) {
		double girth=range.getRange();
		double binsize=girth/numberOfBins;

		int thisBin;
		if (value<=range.getStart()) {
			thisBin=0;
		} else if (value>=range.getStart()+girth) {
			thisBin=numberOfBins-1;
		} else {
			thisBin=(int)Math.ceil((value-range.getStart())/binsize);
		}
		return thisBin;
	}

	public double getMode() {
		double maxProb=-Float.MAX_VALUE;
		int bestIndex=-1;
		for (int i=0; i<histogram.length; i++) {
			if (maxProb<histogram[i]) {
				bestIndex=i;
				maxProb=histogram[i];
			}
		}
		if (bestIndex==-1) {
			return range.getMiddle();
		}
		return bestIndex/(double)numberOfBins*range.getRange()+range.getStart();
	}
	
	@Override
	public double getProbability(double x) {
		return getPDF(x)*getPrior();
	}

	@Override
	public double getPDF(double value) {
		if (sumPriors==0.0f) return 0.0; // no probability

		if (spline!=null) {
			if (spline.isValidPoint(value)) {
				return spline.value(value)/sumPriors;
			} else {
				return 0.0;
			}
		} else {
			double sum=0.0;
			for (Distribution dist : data) {
				double localMin=(double)(dist.getMean()-2.0f*dist.getStdev());
				double localMax=(double)(dist.getMean()+2.0f*dist.getStdev());
				if (value<localMin||value>localMax) {
					continue;
				}

				double probability=(double)dist.getProbability(value);
				if (!Double.isNaN(probability)&&!Double.isInfinite(probability)) {
					sum+=probability;
				}
			}
			return sum/sumPriors;
		}
	}

	@Override
	public double getCDF(double value) {
		if (sumHistogram==0.0f) return 0.0; // no probability

		int bin=getBin(value);
		double sum=0.0;
		for (int i=bin; i<histogram.length; i++) {
			sum+=histogram[i];
		}
		return sum/sumHistogram;
	}

	@Override
	public double getPrior() {
		return prior;
	}

	@Override
	public double getStdev() {
		return stdev;
	}

	@Override
	public double getMean() {
		return mean;
	}
	
	@Override
	public Distribution clone(double mean, double stdev, double prior) {
		throw new EncyclopediaException("Sorry, cannot create a new KDE from just a mean and standard deviation!");
	}
}
