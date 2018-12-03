package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class Histogram implements Distribution {
	private final float start;
	private final float interval;
	private final float[] histogram;
	private final float prior;
	private final float mean;
	private final float stdev;
	
	public static Histogram generateHistogram(float[] x, float[] counts, float interval) {
		assert(x.length==counts.length);
		float min=Float.MAX_VALUE;
		float max=-Float.MAX_VALUE;
		
		for (int i=0; i<x.length; i++) {
			if (min>x[i]) min=x[i];
			if (max<x[i]) max=x[i];
		}
		
		float[] histogram=new float[(int)Math.ceil((max-min)/interval)];
		for (int i=0; i<x.length; i++) {
			int bin=(int)Math.round((x[i]-min)/interval);
			if (bin<0) {
				bin=0;
			} else if (bin>=histogram.length) {
				bin=histogram.length-1;
			}
			histogram[bin]+=counts[i];
		}
		
		return new Histogram(min, interval, histogram);
	}
	
	public float[] getPercentiles(float[] targetPercentiles, float start, float stop) {
		int startBin=getBin(start);
		int stopBin=getBin(stop);
		
		float total=General.sum(histogram, new IntRange(startBin, stopBin));
		float runningSum=0.0f;
		float[] values=new float[targetPercentiles.length];
		Arrays.fill(values, -1.0f);
		values[0]=start;
		values[values.length-1]=stop;
		for (int i=startBin; i<=stopBin; i++) {
			runningSum+=histogram[i]/total;
			for (int j=0; j<values.length; j++) {
				if (values[j]<0&&targetPercentiles[j]<runningSum) {
					values[j]=getValue(i);
				}
			}
		}
		return values;
	}

	public Histogram(float start, float interval, float[] histogram) {
		this(start, interval, histogram, 1.0f);
	}
	public Histogram(float start, float interval, float[] histogram, float prior) {
		this.start=start;
		this.interval=interval;
		this.prior=prior;
		
		float total=General.sum(histogram);
		this.histogram=General.divide(histogram, total);
		
		float totalValue=0.0f;
		for (int i=0; i<this.histogram.length; i++) {
			float y=getValue(i+0.5f);
			totalValue+=y*this.histogram[i];
		}
		mean=totalValue; // this.histogram integration is 1

		float sumSquares=0.0f;
		for (int i=0; i<this.histogram.length; i++) {
			float y=getValue(i+0.5f);
			float delta=y-mean;
			sumSquares+=this.histogram[i]*delta*delta;
		}
		stdev=sumSquares; // this.histogram integration is 1
	}

	@Override
	public String getName() {
		return "histogram";
	}

	@Override
	public Distribution clone(double mean, double stdev, double prior) {
		throw new EncyclopediaException("Sorry, you can't clone a new distribution from a histogram!");
	}

	@Override
	public double getProbability(double x) {
		return getPDF(x)*getPrior();
	}

	@Override
	public double getPDF(double value) {
		int bin=getBin(value);
		return histogram[bin];
	}

	@Override
	public double getCDF(double value) {
		int bin=getBin(value);
		double sum=0.0;
		for (int i=bin; i<histogram.length; i++) {
			sum+=histogram[i];
		}
		return sum;
	}
	
	public int getBin(double value) {
		int bin=(int)Math.round((value-start)/interval);
		if (bin<0) {
			bin=0;
		} else if (bin>=histogram.length) {
			bin=histogram.length-1;
		}
		return bin;
	}
	
	public float getValue(float bin) {
		return start+interval*bin;
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
}
