package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class LinearRegression implements Function {
	public final float m;
	public final float b;

	public final Range xRange;
	public final Range yRange;
	
	public LinearRegression(ArrayList<XYPoint> values) {
		Pair<double[], double[]> xys=XYTrace.toArrays(values);
		float[] x=General.toFloatArray(xys.x);
		float[] y=General.toFloatArray(xys.y);
		
		Pair<Float, Float> regression=getRegression(x, y);
		m=regression.x;
		b=regression.y;
		xRange=new Range(General.min(x), General.max(x));
		yRange=new Range(General.min(y), General.max(y));
	}

	public LinearRegression(float[] x, float[] y) {
		Pair<Float, Float> regression=getRegression(x, y);
		m=regression.x;
		b=regression.y;
		xRange=new Range(General.min(x), General.max(x));
		yRange=new Range(General.min(y), General.max(y));
	}

	@Override
	public float getYValue(float xi) {
		return xi*m+b;
	}

	@Override
	public boolean isXInsideBoundaries(float xi) {
		return xRange.contains(xi);
	}

	@Override
	public boolean isYInsideBoundaries(float yi) {
		return yRange.contains(yi);
	}

	@Override
	public float getXValue(float yi) {
		return (yi-b)/m;
	}

	@Override
	public ArrayList<XYPoint> getKnots() {
		ArrayList<XYPoint> points=new ArrayList<>();
		points.add(new XYPoint(xRange.getStart(), getYValue(xRange.getStart())));
		points.add(new XYPoint(xRange.getStop(), getYValue(xRange.getStop())));
		return points;
	}

	public static Pair<Float, Float> getRegression(float[] x, float[] y) {
		float sumX=0.0f;
		float sumY=0.0f;
		float sumXY=0.0f;
		float sumXX=0.0f;

		for (int i=0; i<y.length; i++) {
			sumX+=x[i];
			sumY+=y[i];
			sumXY+=x[i]*y[i];
			sumXX+=x[i]*x[i];
		}

		float m=((x.length*sumXY)-(sumX*sumY))/((x.length*sumXX)-(sumX*sumX));
		float b=(sumY-m*sumX)/x.length;
		return new Pair<Float, Float>(m, b);
	}
}
