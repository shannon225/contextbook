package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import java.util.ArrayList;

public class WeightedValue extends XYPoint {
	public WeightedValue(double x, double y) {
		super(x, y);
	}

	public double getValue() {
		return getX();
	}
	
	public double getWeight() {
		return getY();
	}
	
	public double average(ArrayList<WeightedValue> values) {
		if (values==null||values.size()==0) return 0.0;
		
		double weightSum=0.0;
		double valueSum=0.0;
		for (WeightedValue v : values) {
			weightSum+=v.getWeight();
			valueSum+=v.getValue()*v.getWeight();
		}
		return valueSum/weightSum;
	}

	public double stdev(ArrayList<WeightedValue> values) {
		if (values==null||values.size()==0) return 0.0;
		
		double weightSum=0.0;
		double valueSum=0.0;
		int nonZeroWeights=0;
		for (WeightedValue v : values) {
			if (v.getWeight()>0.0) {
				nonZeroWeights++;
				weightSum+=v.getWeight();
				valueSum+=v.getValue()*v.getWeight();
			}
		}
		double weightedAverage=valueSum/weightSum;

		double residualSum=0.0;
		for (WeightedValue v : values) {
			double delta=v.getValue()-weightedAverage;
			residualSum+=v.getWeight()*delta*delta;
		}

		double denom=((nonZeroWeights-1)*weightSum)/nonZeroWeights;

		return Math.sqrt(residualSum/denom);
	}
}
