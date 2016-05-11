package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class LinearInterpolatedFunction implements Function {
	private final ArrayList<XYPoint> knots;
	private final double[] x; 
	private final double[] y; 
	public LinearInterpolatedFunction(ArrayList<XYPoint> knots) {
		this.knots=knots;
		Pair<double[], double[]> xys=XYTrace.toArrays(knots);
		this.x=xys.x;
		this.y=xys.y;
	}

	@Override
	public ArrayList<XYPoint> getKnots() {
		return knots;
	}
	
	@Override
	public boolean isXInsideBoundaries(float xi) {
		int upperBin=calculateBinNumber(xi, x);
		if (upperBin==0) return false;
		if (upperBin==x.length) return false;
		return true;
	}

	@Override
	public float getYValue(float xi) {
		int upperBin=calculateBinNumber(xi, x);

		// boundary conditions
		if (upperBin<=0) return (float)y[0];
		if (upperBin>=x.length) return (float)y[y.length-1];

		return linearInterp(x[upperBin-1], (float)xi, x[upperBin], y[upperBin-1], y[upperBin]);
	}
	
	@Override
	public boolean isYInsideBoundaries(float yi) {
		int upperBin=calculateBinNumber(yi, y);
		if (upperBin==0) return false;
		if (upperBin==y.length) return false;
		return true;
	}

	@Override
	public float getXValue(float yi) {
		int upperBin=calculateBinNumber(yi, y);

		// boundary conditions
		if (upperBin==0) return (float)x[0];
		if (upperBin==y.length) return (float)x[y.length-1];

		return linearInterp(y[upperBin-1], (float)yi, y[upperBin], x[upperBin-1], x[upperBin]);
	}
	
	public static float linearInterp(double minX, double X, double maxX, double minY, double maxY) {
		double deltaX=maxX-minX;
		if (deltaX==0) {
			return (float)(maxY+minY)/2f;
		}
		double deltaY=maxY-minY;
		if (deltaY==0) {
			return (float)maxY;
		}
		float interp=(float)(((maxY-minY)/deltaX)*(X-minX)+minY);
		return interp;
	}

	public static int calculateBinNumber(double x, double[] xs) {
		int binNumber=Arrays.binarySearch(xs, x);
		if (binNumber<0) {
			binNumber=(-(binNumber+1));
		}

		return binNumber;
	}
}
