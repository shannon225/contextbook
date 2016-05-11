package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import gnu.trove.list.array.TFloatArrayList;

public class RunningMedianWarper extends LinearInterpolatedFunction {
	public RunningMedianWarper(ArrayList<XYPoint> points, int order, boolean onlyAscending) {
		super(generateKnots(points, order, onlyAscending));
	}

	public static ArrayList<XYPoint> generateKnots(ArrayList<XYPoint> points, int order, boolean onlyAscending) {
		double minX=Double.MAX_VALUE;
		double minY=Double.MAX_VALUE;
		double maxX=-Double.MAX_VALUE;
		double maxY=-Double.MAX_VALUE;
		for (XYPoint xyPoint : points) {
			if (xyPoint.x>maxX) maxX=xyPoint.x;
			if (xyPoint.y>maxY) maxY=xyPoint.y;
			if (xyPoint.x<minX) minX=xyPoint.x;
			if (xyPoint.y<minY) minY=xyPoint.y;
		}
		ArrayList<XYPoint> knots=warp(points, order, onlyAscending);
		if (knots.size()==0) {
			knots.add(new XYPoint(0,0));
		} else {
			knots.add(new XYPoint(minX, minY));
			knots.add(new XYPoint(maxX, maxY));
			Collections.sort(knots);
		}
		return knots;
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
