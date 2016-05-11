package edu.washington.gs.maccoss.encyclopedia.utils.math.distributions;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;

public class TwoDimensionalKDE {
	private final int resolution=1000;
	private final float[][] twoDimensionalHistogram;
	
	private final Range xRange;
	private final Range yRange;

	public TwoDimensionalKDE(ArrayList<XYPoint> points) {
		float[][] stamp=getStamp(points.size());
		
		// first calculate ranges
		float minX=Float.MAX_VALUE;
		float minY=Float.MAX_VALUE;
		float maxX=-Float.MAX_VALUE;
		float maxY=-Float.MAX_VALUE;
		for (XYPoint xyPoint : points) {
			if (xyPoint.x>maxX) maxX=(float)xyPoint.x;
			if (xyPoint.y>maxY) maxY=(float)xyPoint.y;
			if (xyPoint.x<minX) minX=(float)xyPoint.x;
			if (xyPoint.y<minY) minY=(float)xyPoint.y;
		}
		xRange=new Range(minX, maxX);
		yRange=new Range(minY, maxY);

		// then calculate 2d histogram and stamp down the sub-distributions
		twoDimensionalHistogram=new float[resolution][];
		for (int i=0; i<twoDimensionalHistogram.length; i++) {
			twoDimensionalHistogram[i]=new float[resolution];
		}
		
		for (XYPoint xyPoint : points) {
			int xIndex=xRange.linearInterp((float)xyPoint.getX(), 0, resolution-1);
			int yIndex=yRange.linearInterp((float)xyPoint.getY(), 0, resolution-1);
			stampDistribution(xIndex, yIndex, stamp);
		}

	}
	
	public Function trace() {
		float max=0f;
		int maxXIndex=0;
		int maxYIndex=0;
		for (int i=0; i<twoDimensionalHistogram.length; i++) {
			for (int j=0; j<twoDimensionalHistogram[i].length; j++) {
				if (twoDimensionalHistogram[i][j]>max) {
					max=twoDimensionalHistogram[i][j];
					maxXIndex=i;
					maxYIndex=j;
				}
			}
		}
		
		ArrayList<XYPoint> points=new ArrayList<XYPoint>();
		traceSouthWest(maxXIndex, maxYIndex, points);
		traceNorthEast(maxXIndex, maxYIndex, points);
		Collections.sort(points);
		
		return new LinearInterpolatedFunction(points);
	}
	
	public void traceNorthEast(int i, int j, ArrayList<XYPoint> prev) {
		float x=xRange.mapBackToRange(i, 0, resolution-1);
		float y=yRange.mapBackToRange(j, 0, resolution-1);
		prev.add(new XYPoint(x, y));
		
		if (i>=resolution-1||j>=resolution-1) {
			return;
		}
		
		float east=twoDimensionalHistogram[i+1][j];
		float northeast=twoDimensionalHistogram[i+1][j+1];
		float north=twoDimensionalHistogram[i][j+1];
		float max=Math.max(Math.max(east, north), northeast);
		
		if (northeast==max||east==north) {
			traceNorthEast(i+1, j+1, prev);
		} else if (east==max) {
			traceNorthEast(i+1, j, prev);
		} else {
			traceNorthEast(i, j+1, prev);
		}
	}
	
	public void traceSouthWest(int i, int j, ArrayList<XYPoint> prev) {
		float x=xRange.mapBackToRange(i, 0, resolution-1);
		float y=yRange.mapBackToRange(j, 0, resolution-1);
		prev.add(new XYPoint(x, y));
		
		if (i<=0||j<=0) {
			return;
		}
		
		float west=twoDimensionalHistogram[i-1][j];
		float southwest=twoDimensionalHistogram[i-1][j-1];
		float south=twoDimensionalHistogram[i][j-1];
		float max=Math.max(Math.max(west, south), southwest);
		
		if (southwest==max||west==south) {
			traceSouthWest(i-1, j-1, prev);
		} else if (west==max) {
			traceSouthWest(i-1, j, prev);
		} else {
			traceSouthWest(i, j-1, prev);
		}
	}
	
	public void plot() {
		ArrayList<XYZPoint> heatData=new ArrayList<XYZPoint>();
		for (int i=0; i<twoDimensionalHistogram.length; i++) {
			float x=xRange.mapBackToRange(i, 0, resolution-1);
			for (int j=0; j<twoDimensionalHistogram[i].length; j++) {
				float y=yRange.mapBackToRange(j, 0, resolution-1);
				heatData.add(new XYZPoint(x, y, twoDimensionalHistogram[i][j]));
			}
		}
		Charter.launchChart("X", "Y", false, new XYZTrace("Density", heatData));
	}
	
	public float[][] getStamp(int n) {
		Distribution dist=new CosineGaussian(0.0f, Math.pow(n, 1f/3f), 1.0f);
		int stampRadius=Math.round(2.0f*(float)dist.getStdev());
		
		float[][] stamp=new float[stampRadius*2+1][];
		for (int i=0; i<stamp.length; i++) {
			stamp[i]=new float[stampRadius*2+1];
			for (int j=0; j<stamp[i].length; j++) {
				int deltaX=Math.abs(i-stampRadius);
				int deltaY=Math.abs(j-stampRadius);
				double distance=Math.sqrt(deltaX*deltaX+deltaY*deltaY);
				stamp[i][j]=(float)dist.getPDF(distance);
			}
		}
		return stamp;
	}

	public void stampDistribution(int xIndex, int yIndex, float[][] stamp) {
		int stampRadius=stamp.length/2;
		for (int i=0; i<stamp.length; i++) {
			int localX=xIndex+i-stampRadius;
			if (localX>=0&&localX<twoDimensionalHistogram.length) {
				for (int j=0; j<stamp[i].length; j++) {
					int localY=yIndex+j-stampRadius;
					if (localY>=0&&localY<twoDimensionalHistogram[i].length) {
						twoDimensionalHistogram[localX][localY]+=stamp[i][j];
					}
				}
			}
		}
	}
}
