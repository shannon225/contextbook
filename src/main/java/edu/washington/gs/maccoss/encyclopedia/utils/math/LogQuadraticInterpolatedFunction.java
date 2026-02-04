package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class LogQuadraticInterpolatedFunction extends PeakShapeFunction {

    private final ArrayList<XYPoint> knots;
    private final double[] x;
    private final double[] y;
    private final double averageXIncrement;
    private final double minimumMeaningfulValue;

	public LogQuadraticInterpolatedFunction(ArrayList<XYPoint> knots) {
		this(knots, Math.E);
	}
	public LogQuadraticInterpolatedFunction(ArrayList<XYPoint> knots, double minimumMeaningfulValue) {
		Pair<double[], double[]> xys=XYTrace.toArrays(knots);
		this.minimumMeaningfulValue=minimumMeaningfulValue;
		this.x=xys.x;
		this.y=Log.logLn(General.add(xys.y, this.minimumMeaningfulValue));
		this.knots=XYTrace.toPoints(x, y);
		if (knots.size()>1) {
			this.averageXIncrement=(x[x.length-1]-x[0])/x.length;
		} else {
			this.averageXIncrement=0.0;
		}
	}

	public LogQuadraticInterpolatedFunction(double[] x, double[] y) {
		this(x, y, Math.E);
	}
	public LogQuadraticInterpolatedFunction(double[] x, double[] y, double minimumMeaningfulValue) {
		this.minimumMeaningfulValue=minimumMeaningfulValue;
		this.x=x;
		this.y=Log.logLn(General.add(y, this.minimumMeaningfulValue));
		this.knots=XYTrace.toPoints(x, y);
		this.averageXIncrement=(x[x.length-1]-x[0])/x.length;
	}

    @Override
    public ArrayList<XYPoint> getKnots() {
        return knots;
    }

    @Override
    public boolean isXInsideBoundaries(float xi) {
        int upperBin = calculateBinNumber(xi, x);
        if (upperBin == 0) return false;
        if (upperBin == x.length) return false;
        return true;
    }

    @Override
    public float getYValue(float xi) {
		if (x.length==1) return 0.0f;
		if (xi<=x[0]) return 0.0f;
		if (xi>=x[x.length-1]) return 0.0f;

		// binary search for insertion point
		int lo=0, hi=x.length-1;
		while (lo+1<hi) {
			int mid=lo+(hi-lo)/2;
			if (x[mid]<=xi) {
				lo=mid;
			} else {
				hi=mid;
			}
		}
		int i=lo;
		
		double y0=i>0?y[i-1]:0.0f;
		double y1=i>=0?y[i]:0.0f;
		double y2=i<y.length-1?y[i+1]:0.0f;
		double y3=i<y.length-2?y[i+2]:0.0f;

		double x0=i>0?x[i-1]:x[i]-averageXIncrement;
		double x1=i>=0?x[i]:0.0f;
		double x2=i<x.length-1?x[i+1]:x[i]+averageXIncrement;
		double x3=i<x.length-2?x[i+2]:x[i]+2.0*averageXIncrement;

		if (y1==0.0f&&y2==0.0f) {
			return 0.0f;
		}

		double logIntensity=Math.max(0.0f, quadraticFit(x0, x1, x2, x3, y0, y1, y2, y3, xi));
		float intensity=(float)(Math.exp(logIntensity)-minimumMeaningfulValue);
		if (intensity<0.0f) return 0.0f; 
		return intensity;
    }

    @Override
    public boolean isYInsideBoundaries(float yi) {
        int upperBin = calculateBinNumber(yi, y);
        if (upperBin == 0) return false;
        if (upperBin == y.length) return false;
        return true;
    }

    @Override
    public float getXValue(float yi) {
    	throw new EncyclopediaException("Quadratic cannot be solved for X since it produces zero or two values!");
    }

    /** Binary search boundary index like the linear version. */
    public static int calculateBinNumber(double x, double[] xs) {
        int binNumber = Arrays.binarySearch(xs, x);
        if (binNumber < 0) binNumber = (-(binNumber + 1));
        return binNumber;
	}

	/**
	 * Calculates a quadratic fit forced to go between the 
	 * middle points, ensuring continuity. Basic fallbacks 
	 * exist for if x2-x1=0 and if all the points fall on a 
	 * straight line to avoid micro-oscillations.
	 * 
	 * Calculated as a linear portion plus the parabolic 
	 * basis. This approach implicitly forces the fit   
	 * through the points x1,y1 and x2,y2:
	 * 
	 *    y = slope*xi+b + k*(x-x1)*(x-x2)
	 *    
	 * This calculates a line between the inner points and 
	 * adjusts it with the parabolic shape. K is how strong 
	 * the adjustment factor is.
	 *    
	 * Converting to y=a*x^2+b*x+c
	 *    a = k
	 *    b = slope-k(x1+x2)
	 *    c = b+k*x1*x2
	 */
	static double quadraticFit(
			double x0, double x1, double x2, double x3, // x values surrounding xi
			double y0, double y1, double y2, double y3, // corresponding y values
			double xi) {
		
		// If inner times coincide, fall back to their average
		double epsilon=1e-6;
		double deltaX12=x2-x1;
		if (Math.abs(deltaX12)<epsilon) return 0.5*(y1+y2);

		// Calculate a line through inner points
		double slope=(y2-y1)/deltaX12;
		double b=y1-slope*x1;

		// Check to see if just effectively a straight line within epsilon
		double linearFitAtX0=slope*x0+b;
		double linearFitAtX3=slope*x3+b;
		double linearFitAtXi=slope*xi+b;
		if (Math.abs(linearFitAtX0-y0)<=epsilon&&Math.abs(linearFitAtX3-y3)<=epsilon) {
			// If so, just return the line
			return linearFitAtXi;
		}

		// Calculate the parabolic shape relative to xi
		double parabolaXi=(xi-x1)*(xi-x2);

		// Calculate the parabolic shape at the outer points as left and right guides
		double parabolaX0=(x0-x1)*(x0-x2);
		double parabolaX3=(x3-x1)*(x3-x2);

		// Linear weights that defer to higher outer values
		double weight0=Math.max(y0, epsilon);
		double weight3=Math.max(y3, epsilon);

		// Closed-form k from weighted least squares using the two outer points
		double numerator=weight0*parabolaX0*(y0-linearFitAtX0)+weight3*parabolaX3*(y3-linearFitAtX3);
		double denominator=weight0*parabolaX0*parabolaX0+weight3*parabolaX3*parabolaX3;

		// Final curvature parameter chosen by the weighted fit
		double k=(Math.abs(denominator)>epsilon)?(numerator/denominator):0.0;
		
		//System.out.println(xi+"\t"+k+"\t"+(slope-k*(x1+x2))+"\t"+(b+k*x1*x2));

		// Quadratic value at rtInSec
		return linearFitAtXi+k*parabolaXi;
	}

}