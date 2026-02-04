package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class LinearInterpolatedFunction extends PeakShapeFunction {
	private final ArrayList<XYPoint> knots;
	private final double[] x; 
	private final double[] y; 
	
	public LinearInterpolatedFunction(ArrayList<XYPoint> knots) {
		this.knots=knots;
		Pair<double[], double[]> xys=XYTrace.toArrays(knots);
		this.x=xys.x;
		this.y=xys.y;
	}
	
	public LinearInterpolatedFunction(double[] x, double[] y) {
		this.x=x;
		this.y=y;
		knots=XYTrace.toPoints(x, y);
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

	@Override
	public XYPoint getApex(double left, double right) {
		int leftBin=calculateBinNumber(left, x);
		int rightBin=calculateBinNumber(right, x);
		
		int maxBin=leftBin;
		for (int i=leftBin+1; i<=Math.min(y.length-1, rightBin); i++) {
			if (y[i]>y[maxBin]) {
				maxBin=i;
			}
		}
		return new XYPoint(x[maxBin], y[maxBin]);
	}

	@Override
	public double integrate(double left, double right) {
	    if (x == null || y == null || x.length < 2) return 0.0;
	    if (right <= left) return 0.0;

	    // Clip to data domain [x[0], x[n-1]]
	    final double a = Math.max(left, x[0]);
	    final double b = Math.min(right, x[x.length - 1]);
	    if (b <= a) return 0.0;

	    // Find bins that contain the clipped endpoints
	    int binStart = calculateBinNumber(a, x);   // index of right endpoint of the segment containing 'a'
	    int binEnd   = calculateBinNumber(b, x);   // index of right endpoint of the segment containing 'b'

	    // Ensure we loop over valid segments [i-1, i]
	    int iStart = Math.max(1, binStart);
	    int iEnd   = Math.min(binEnd, x.length - 1);

	    double total = 0.0;

	    for (int i = iStart; i <= iEnd; i++) {
	        final double xl = Math.max(a, x[i - 1]);
	        final double xr = Math.min(b, x[i]);
	        if (xr <= xl) continue;

	        // Linear values at xl and xr within segment [x[i-1], x[i]]
	        final double xLseg0 = x[i - 1];
	        final double xRseg0 = x[i];
	        final double yLseg0 = y[i - 1];
	        final double yRseg0 = y[i];

	        double yL, yR;
	        if (xRseg0 == xLseg0) {
	            // Degenerate segment, treat as constant
	            yL = 0.5 * (yLseg0 + yRseg0);
	            yR = yL;
	        } else {
	            double invDx = 1.0 / (xRseg0 - xLseg0);
	            // y(x) = yL0 + (yR0 - yL0) * (x - xL0)/(xR0 - xL0)
	            double slope = (yRseg0 - yLseg0) * invDx;
	            yL = yLseg0 + slope * (xl - xLseg0);
	            yR = yLseg0 + slope * (xr - xLseg0);
	        }

	        total += 0.5 * (yL + yR) * (xr - xl);
	    }

	    return total;
	}

	private static int calculateBinNumber(double x, double[] xs) {
		int binNumber=Arrays.binarySearch(xs, x);
		if (binNumber<0) {
			binNumber=(-(binNumber+1));
		}

		return binNumber;
	}
}
