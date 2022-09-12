package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.CosineGaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import gnu.trove.list.array.TFloatArrayList;

public class RectangleKDE {
	private static final double BANDWIDTH_TO_STDEV = 2.0*Math.sqrt(2.0*Math.log(2.0));
	private final int yResolution;
	private final int xResolution;
	private final float[][] twoDimensionalHistogram;
	
	private final Range xRange;
	private final Range yRange;
	
	public RectangleKDE(ArrayList<XYPoint> points, int xResolution, int yResolution) {
		this.yResolution=yResolution;
		this.xResolution=xResolution;
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

		TFloatArrayList xs=new TFloatArrayList();
		TFloatArrayList ys=new TFloatArrayList();
		for (XYPoint xyPoint : points) {
			int xIndex=xRange.linearInterp((float)xyPoint.getX(), 0, xResolution-1);
			int yIndex=yRange.linearInterp((float)xyPoint.getY(), 0, yResolution-1);
			xs.add(xIndex);
			ys.add(yIndex);
		}
		
		float xstdev=General.stdev(xs.toArray());
		float ystdev=General.stdev(ys.toArray());
		
		// Silverman's rule of thumb
		double bandwidthX=Math.pow(points.size(), -1f/6f)*xstdev;
		// division by 2.3548 converts bandwidth (fwhm) to stdev for gaussians
		double stdevX=Math.min(xResolution/10d/4d, bandwidthX/BANDWIDTH_TO_STDEV); // max of 10% distributions (4x stdev)
		Distribution distX=new CosineGaussian(0.0, stdevX, 1.0);
		
		// Silverman's rule of thumb
		double bandwidthY=Math.pow(points.size(), -1f/6f)*ystdev;
		// division by 2.3548 converts bandwidth (fwhm) to stdev for gaussians
		double stdevY=Math.min(yResolution/10d/4d, bandwidthY/BANDWIDTH_TO_STDEV); // max of 10% distributions (4x stdev)
		Distribution distY=new CosineGaussian(0.0, stdevY, 1.0);
		
		float[][] stamp=getStamp(distX, distY);
		
		// then calculate 2d histogram and stamp down the sub-distributions
		twoDimensionalHistogram=new float[xResolution][];
		for (int i=0; i<twoDimensionalHistogram.length; i++) {
			twoDimensionalHistogram[i]=new float[yResolution];
		}
		
		for (XYPoint xyPoint : points) {
			int xIndex=xRange.linearInterp((float)xyPoint.getX(), 0, xResolution-1);
			int yIndex=yRange.linearInterp((float)xyPoint.getY(), 0, yResolution-1);
			stampDistribution(xIndex, yIndex, stamp);
		}

	}

	public float[][] getTwoDimensionalHistogram() {
		return twoDimensionalHistogram;
	}
	public Range getXRange() {
		return xRange;
	}
	
	public Range getYRange() {
		return yRange;
	}
	
	public double f(double X, double Y) {
		int xIndex=xRange.linearInterp((float)X, 0, xResolution-1);
		int yIndex=yRange.linearInterp((float)Y, 0, yResolution-1);
		if (xIndex<0||xIndex>=xResolution) return 0.0;
		if (yIndex<0||yIndex>=yResolution) return 0.0;
		
		return twoDimensionalHistogram[xIndex][yIndex];
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
		traceWest(maxXIndex, maxYIndex, points);
		traceEast(maxXIndex, maxYIndex, points);
		Collections.sort(points);
		
		return new LinearInterpolatedFunction(points);
	}
	
	public void traceEast(int i, int j, ArrayList<XYPoint> prev) {
		float x=xRange.mapBackToRange(i, 0, xResolution-1);
		float y=yRange.mapBackToRange(j, 0, yResolution-1);
		prev.add(new XYPoint(x, y));
		
		if (i>=xResolution-1) return;
		if (j<=0) {
			j=1;
		} else if (j>=yResolution-1) {
			j=yResolution-2;
		}
		
		float east=twoDimensionalHistogram[i+1][j];
		float northeast=twoDimensionalHistogram[i+1][j+1];
		float southeast=twoDimensionalHistogram[i+1][j-1];
		float max=Math.max(Math.max(east, southeast), northeast);
		
		if (east==max||northeast==southeast) {
			traceEast(i+1, j, prev);
		} else if (northeast==max) {
			traceEast(i+1, j+1, prev);
		} else {
			traceEast(i+1, j-1, prev);
		}
	}
	
	public void traceWest(int i, int j, ArrayList<XYPoint> prev) {
		float x=xRange.mapBackToRange(i, 0, xResolution-1);
		float y=yRange.mapBackToRange(j, 0, yResolution-1);
		prev.add(new XYPoint(x, y));
		
		if (i<=0) return;
		if (j<=0) {
			j=1;
		} else if (j>=yResolution-1) {
			j=yResolution-2;
		}
		
		float west=twoDimensionalHistogram[i-1][j];
		float southwest=twoDimensionalHistogram[i-1][j-1];
		float northwest=twoDimensionalHistogram[i-1][j+1];
		float max=Math.max(Math.max(west, northwest), southwest);
		
		if (west==max||southwest==northwest) {
			traceWest(i-1, j, prev);
		} else if (southwest==max) {
			traceWest(i-1, j-1, prev);
		} else {
			traceWest(i-1, j+1, prev);
		}
	}
	
	public void plot() {
		ArrayList<XYZPoint> heatData=new ArrayList<XYZPoint>();
		for (int i=0; i<twoDimensionalHistogram.length; i++) {
			float x=xRange.mapBackToRange(i, 0, xResolution-1);
			for (int j=0; j<twoDimensionalHistogram[i].length; j++) {
				float y=yRange.mapBackToRange(j, 0, yResolution-1);
				heatData.add(new XYZPoint(x, y, twoDimensionalHistogram[i][j]));
			}
		}
		Charter.launchChart("X", "Y", false, new XYZTrace("Density", heatData));
	}
	
	public static float[][] getStamp(Distribution distX, Distribution distY) {
		int stampRadiusX=Math.round(2.0f*(float)distX.getStdev());
		int stampRadiusY=Math.round(2.0f*(float)distY.getStdev());
		
		float[][] stamp=new float[stampRadiusX*2+1][];
		for (int i=0; i<stamp.length; i++) {
			stamp[i]=new float[stampRadiusY*2+1];
			for (int j=0; j<stamp[i].length; j++) {
				int deltaX=Math.abs(i-stampRadiusX);
				int deltaY=Math.abs(j-stampRadiusY);
				stamp[i][j]=(float)Math.sqrt(distX.getPDF(deltaX)*distY.getPDF(deltaY));
			}
		}
		return stamp;
	}

	public void stampDistribution(int xIndex, int yIndex, float[][] stamp) {
		int stampRadiusX=stamp.length/2;
		for (int i=0; i<stamp.length; i++) {
			int localX=xIndex+i-stampRadiusX;
			int stampRadiusY=stamp[i].length/2;
			if (localX>=0&&localX<twoDimensionalHistogram.length) {
				for (int j=0; j<stamp[i].length; j++) {
					if (stamp[i][j]>0.0f) {
						int localY=yIndex+j-stampRadiusY;
						if (localY>=0&&localY<twoDimensionalHistogram[i].length) {
							twoDimensionalHistogram[localX][localY]+=stamp[i][j];
						}
					}
				}
			}
		}
	}
}
