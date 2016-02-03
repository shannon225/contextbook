package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import gnu.trove.list.array.TFloatArrayList;

public class RunningMedianWarper implements Function {
	private final ArrayList<XYPoint> knots;
	private final double[] x; 
	private final double[] y; 
	public RunningMedianWarper(ArrayList<XYPoint> points, int order, boolean onlyAscending) {
		knots=warp(points, order, onlyAscending);
		if (knots.size()==0) {
			knots.add(new XYPoint(0,0));
			this.x=new double[1];
			this.y=new double[1];
		} else {
			Collections.sort(knots);
			Pair<double[], double[]> xys=XYTrace.toArrays(knots);
			this.x=xys.x;
			this.y=xys.y;
		}
	}
	
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
	
	public static ArrayList<XYPoint> warp(ArrayList<XYPoint> points, int order, boolean onlyAscending) {
		if (points.size()==0) return points;
		
		Pair<float[], float[]>xys=XYTrace.toFloatArrays(points);
		float[] x=xys.x;
		float[] y=xys.y;
		float[] newX=new float[x.length];
		float[] newY=new float[x.length];
		for (int i=0; i<y.length; i++) {
			newX[i]=getMedian(x, order, i);
			newY[i]=getMedian(y, order, i);
		}
		ArrayList<XYPoint> values=new ArrayList<XYPoint>();
		float lastX=newX[newY.length/2];
		float lastY=newY[newY.length/2];
		for (int i=newY.length/2; i<newY.length; i++) {
			if (onlyAscending) {
				if (newX[i]<lastX) continue;
				if (newY[i]<lastY) continue;
			}
			values.add(new XYPoint(newX[i], newY[i]));
			lastX=newX[i];
			lastY=newY[i];
		}
		
		lastX=newX[newY.length/2];
		lastY=newY[newY.length/2];
		for (int i=newY.length/2-1; i>=0; i--) {
			if (onlyAscending) {
				if (newX[i]>lastX) continue;
				if (newY[i]>lastY) continue;
			}
			values.add(0, new XYPoint(newX[i], newY[i]));
			lastX=newX[i];
			lastY=newY[i];
		}
		
		return values;
	}

	private static float getMedian(float[] y, int order, int i) {
		float[] array=extractRange(y, order, i);
		float median=QuickMedian.median(array);
		
		return median;
	}

	public static float[] extractRange(float[] y, int order, int i) {
		TFloatArrayList selected=new TFloatArrayList();
		int minRange=i-Math.max(0, i-order);
		int maxRange=Math.min(y.length, i+order)-i;
		
		int finalRange=Math.min(minRange, maxRange);
		
		for (int j=i; j>i-finalRange; j--) {
			selected.add(y[j]);
		}
		for (int j=i+1; j<i+finalRange; j++) {
			selected.add(y[j]);
		}
		float[] array=selected.toArray();
		return array;
	}

}
