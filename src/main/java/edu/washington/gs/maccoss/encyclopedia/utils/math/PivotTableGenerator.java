package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public class PivotTableGenerator {
	public static ArrayList<XYPoint> createPivotTable(float[] data) {
		float actualMin=Float.MAX_VALUE;
		float actualMax=-Float.MAX_VALUE;
		
		for (int i=0; i<data.length; i++) {
			if (data[i]>actualMax) {
				actualMax=data[i];
			}
			if (data[i]<actualMin) {
				actualMin=data[i];
			}
		}
		
		int binCount=data.length/20;
		
		return createPivotTable(data, actualMin, actualMax, (actualMax-actualMin)/binCount);
	}
	public static ArrayList<XYPoint> createPivotTable(float[] data, float minimum, float maximum, float binsize) {
		int numberOfBins=(int)((maximum-minimum)/binsize)+1;
		int[] histogram=new int[numberOfBins];
		for (int i=0; i<data.length; i++) {
			int index=Math.round((data[i]-minimum)/binsize);
			if (index<0) index=0;
			if (index>=numberOfBins) index=numberOfBins-1;
			
			histogram[index]++;
		}
		
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		for (int i=0; i<histogram.length; i++) {
			float binCenter=minimum+binsize*(i+0.5f);
			points.add(new XYPoint(binCenter, histogram[i]));
		}
		return points;
	}

}
