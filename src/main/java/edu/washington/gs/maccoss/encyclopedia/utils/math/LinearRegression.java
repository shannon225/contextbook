package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public class LinearRegression {
	public static Pair<Float, Float> getRegression(float[] x, float[] y) {
		float sumX=0.0f;
		float sumY=0.0f;
		float sumXY=0.0f;
		float sumXX=0.0f;
		
		for (int i = 0; i < y.length; i++) {
			sumX+=x[i];
			sumY+=y[i];
			sumXY+=x[i]*y[i];
			sumXX+=x[i]*x[i];
		}
		
		float m=((x.length*sumXY)-(sumX*sumY)) / ((x.length*sumXX)-(sumX*sumX));
		float b=(sumY-m*sumX)/x.length;
		return new Pair<Float, Float>(m, b);
	}
}
