package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.WeightedValue;

public class KDE {
	private final ArrayList<Distribution> data;
	private final Range range;
	private final int numberOfBins=100;
	private final float[] histogram;
	private final float sumPriors;

	public KDE(ArrayList<WeightedValue> values) {
		double stdev=WeightedValue.stdev(values);
		// Silverman's (1986) rule of thumb (wikipedia)
		double bandwidth=stdev*Math.pow(4.0/3.0/values.size(), 1.0/5.0);
		
		data=new ArrayList<Distribution>();
		for (WeightedValue value : values) {
			data.add(new CosineGaussian(value.getValue(), bandwidth, value.getWeight()));
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
		float binsize=(max-min)/numberOfBins;

		histogram=new float[numberOfBins];
		float[] binValues=new float[histogram.length];
		for (int i=0; i<histogram.length; i++) {
			binValues[i]=i*binsize+min;
		}

		float total=0.0f;
		for (Distribution dist : data) {
			total+=dist.getPrior();
			float localMin=(float)(dist.getMean()-2.0f*dist.getStdev());
			float localMax=(float)(dist.getMean()+2.0f*dist.getStdev());
			int startIndex=Math.max(0, (int)Math.floor((localMin-min)/binsize));
			int stopIndex=Math.min(histogram.length-1, (int)Math.ceil((localMax-min)/binsize));

			for (int i=startIndex; i<=stopIndex; i++) {
				float probability=(float)dist.getProbability(binValues[i]);
				if (!Float.isNaN(probability)&&!Float.isInfinite(probability)) {
					histogram[i]+=probability;
				}
			}
		}
		sumPriors=total;
	}

	public int getBin(float value) {
		float girth=range.getRange();
		float binsize=girth/numberOfBins;

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
	
	public float getMode() {
		float maxProb=-Float.MAX_VALUE;
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
		return bestIndex/(float)numberOfBins*range.getRange()+range.getStart();
	}

	public float getProbability(float value) {
		if (sumPriors==0.0f) return 0.0f; // no probability
		
		float sum=0.0f;
		for (Distribution dist : data) {
			float localMin=(float)(dist.getMean()-2.0f*dist.getStdev());
			float localMax=(float)(dist.getMean()+2.0f*dist.getStdev());
			if (value>localMin||value<localMax) {
				continue;
			}
			
			float probability=(float)dist.getProbability(value);
			if (!Float.isNaN(probability)&&!Float.isInfinite(probability)) {
				sum+=probability;
			}
		}
		return sum/sumPriors;
	}
}
