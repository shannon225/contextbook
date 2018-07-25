package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import gnu.trove.set.hash.TDoubleHashSet;

public class PivotTableGenerator {
	public static ArrayList<XYPoint>[] createPivotTables(float[][] datas, boolean removeNonZero) {
		float actualMin=Float.MAX_VALUE;
		float actualMax=-Float.MAX_VALUE;
		
		int length=0;
		for (int a=0; a<datas.length; a++) {
			float[] data=datas[a];
			length+=data.length;
			for (int i=0; i<data.length; i++) {
				if (data[i]>actualMax) {
					actualMax=data[i];
				}
				if (data[i]<actualMin) {
					actualMin=data[i];
				}
			}
		}
		
		int binCount=(length/datas.length)/200;
		binCount=Math.max(binCount, 50);

		@SuppressWarnings("unchecked")
		ArrayList<XYPoint>[] traces=new ArrayList[datas.length];
		for (int a=0; a<datas.length; a++) {
			float[] data=datas[a];
			ArrayList<XYPoint> trace=createPivotTable(data, actualMin, actualMax, (actualMax-actualMin)/binCount);
			traces[a]=trace;
		}

		if (removeNonZero) {
			TDoubleHashSet xs=new TDoubleHashSet();
			for (int i=0; i<traces.length; i++) {
				ArrayList<XYPoint> trace=traces[i];
				for (XYPoint point : trace) {
					if (point.getY()>0.0) {
						xs.add(point.getX());
					}
				}
			}
			for (int i=0; i<traces.length; i++) {
				ArrayList<XYPoint> trace=traces[i];
				ArrayList<XYPoint> keepers=new ArrayList<>();
				for (XYPoint point : trace) {
					if (xs.contains(point.x)) {
						keepers.add(point);
					}
				}
				traces[i]=keepers;
			}
		}
		
		return traces;
	}
	
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
		
		int binCount=data.length/25;
		
		return createPivotTable(data, actualMin, actualMax, (actualMax-actualMin)/binCount);
	}
	public static ArrayList<XYPoint> createPivotTable(float[] data, float minimum, float maximum, float binsize) {
		int numberOfBins=(int)((maximum-minimum)/binsize)+1;
		int[] histogram=new int[numberOfBins];
		for (int i=0; i<data.length; i++) {
			int index=(int)Math.floor((data[i]-minimum)/binsize);
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
