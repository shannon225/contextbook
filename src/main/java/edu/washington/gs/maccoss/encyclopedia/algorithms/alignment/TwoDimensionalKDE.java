package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.gui.general.Charter;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserConsumer;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserProducer;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearInterpolatedFunction;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.CosineGaussian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import gnu.trove.list.array.TFloatArrayList;

public class TwoDimensionalKDE {
	public static final int DEFAULT_RESOLUTION = 1000;
	public static final int HIGHER_RESOLUTION = 3000;
	private static final double BANDWIDTH_TO_STDEV = 2.0*Math.sqrt(2.0*Math.log(2.0));
	private final int resolution;
	private final float[][] twoDimensionalHistogram;
	
	private final Range xRange;
	private final Range yRange;
	
	public static void main(String[] args) {
		File f=new File("/Users/searleb/Documents/projects/encyclopedia/retention_times.txt");

		final ArrayList<XYPoint> rts=new ArrayList<XYPoint>();

		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				float predicted=Float.parseFloat(row.get("predicted"));
				float actual=Float.parseFloat(row.get("actual"));
				rts.add(new XYPoint(predicted, actual));
			}
			
			@Override
			public void cleanup() {
			}
		};

		BlockingQueue<Map<String, String>> blockingQueue=new LinkedBlockingQueue<Map<String, String>>();
		TableParserProducer producer=new TableParserProducer(blockingQueue, f, "\t", 1);
		TableParserConsumer consumer=new TableParserConsumer(blockingQueue, muscle);

		Thread producerThread=new Thread(producer);
		Thread consumerThread=new Thread(consumer);
		producerThread.start();
		consumerThread.start();

		try {
			producerThread.join();
			consumerThread.join();
		} catch (InterruptedException ie) {
			Logger.errorLine("Percolator reading interrupted!");
			Logger.errorException(ie);
		}

		TwoDimensionalKDE filter=null;
		for (int i=0; i<20; i++) {
			long time=System.currentTimeMillis();
			filter=new TwoDimensionalKDE(rts);
			filter.trace();
			System.out.println((System.currentTimeMillis()-time));
		}

//		final TwoDimensionalKDE finalFilter=filter;
//		Mapper mapper=new Mapper() {
//			@Override
//			public double f(double arg0, double arg1) {
//				return finalFilter.f(arg0, arg1);
//			}
//		};
//		Charter3d.plot(mapper, 
//				new org.jzy3d.maths.Range(finalFilter.getXRange().getStart(), finalFilter.getXRange().getStop()), 
//				new org.jzy3d.maths.Range(finalFilter.getYRange().getStart(), finalFilter.getYRange().getStop()), 
//				finalFilter.getResolution()/5);
	}

	public TwoDimensionalKDE(ArrayList<XYPoint> points) {
		this(points, DEFAULT_RESOLUTION);
	}
	
	public TwoDimensionalKDE(ArrayList<XYPoint> points, int resolution) {
		this.resolution=resolution;
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
			int xIndex=xRange.linearInterp((float)xyPoint.getX(), 0, resolution-1);
			int yIndex=yRange.linearInterp((float)xyPoint.getY(), 0, resolution-1);
			xs.add(xIndex);
			ys.add(yIndex);
		}
		
		float xstdev=General.stdev(xs.toArray());
		float ystdev=General.stdev(ys.toArray());
		
		// Silverman's rule of thumb
		double bandwidth=Math.pow(points.size(), -1f/6f)*(xstdev+ystdev)/2.0f;
		// division by 2.3548 converts bandwidth (fwhm) to stdev for gaussians
		double stdev=Math.min(resolution/10d/4d, bandwidth/BANDWIDTH_TO_STDEV); // max of 10% distributions (4x stdev)
		Distribution dist=new CosineGaussian(0.0, stdev, 1.0);
		float[][] stamp=getStamp(dist);
		
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

	public float[][] getTwoDimensionalHistogram() {
		return twoDimensionalHistogram;
	}
	public Range getXRange() {
		return xRange;
	}
	
	public Range getYRange() {
		return yRange;
	}
	
	public int getResolution() {
		return resolution;
	}
	
	public double f(double X, double Y) {
		int xIndex=xRange.linearInterp((float)X, 0, resolution-1);
		int yIndex=yRange.linearInterp((float)Y, 0, resolution-1);
		if (xIndex<0||xIndex>=resolution) return 0.0;
		if (yIndex<0||yIndex>=resolution) return 0.0;
		
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
		traceSouthWest(maxXIndex, maxYIndex, points);
		traceNorthEast(maxXIndex, maxYIndex, points);
		Collections.sort(points);
		
		return new LinearInterpolatedFunction(points);
	}
	
	public void traceNorthEast(int initialI, int initialJ, ArrayList<XYPoint> prev) {
		int i = initialI;
		int j = initialJ;

		while (true) {// gauranteed to exit because each step either i or j (or both) increases
			float x = xRange.mapBackToRange(i, 0, resolution - 1);
			float y = yRange.mapBackToRange(j, 0, resolution - 1);
			prev.add(new XYPoint(x, y));

			if (i >= resolution - 1 || j >= resolution - 1) {
				return;
			}

			float east = twoDimensionalHistogram[i + 1][j];
			float northeast = twoDimensionalHistogram[i + 1][j + 1];
			float north = twoDimensionalHistogram[i][j + 1];
			float max = Math.max(Math.max(east, north), northeast);

			if (northeast == max || east == north) {
				i++;
				j++;
			} else if (east == max) {
				// We want to go east
				// First check if the last move was north
				// (but only if there was a previous move)
				if (prev.size() >= 2) {
					final XYPoint last = prev.get(prev.size() - 2);
					if (last.x == x) {
						// If the previous points X coord is the same
						// as this point, the last move was to the north.
						// Thus we remove this current point and just move
						// northeast from the last point.
						prev.remove(prev.size() - 1);
					}
				}
				i++;
			} else {
				// We want to go north
				// First check if the last move was east
				// (but only if there was a previous move)
				if (prev.size() >= 2) {
					final XYPoint last = prev.get(prev.size() - 2);
					if (last.y == y) {
						// If the previous points Y coord is the same
						// as this point, the last move was to the east.
						// Thus we remove this current point and just move
						// northeast from the last point.
						prev.remove(prev.size() - 1);
					}
				}
				j++;
			}
		}
	}
	
	public void traceSouthWest(int initialI, int initialJ, ArrayList<XYPoint> prev) {
		int i = initialI;
		int j = initialJ;

		while (true) { // gauranteed to exit because each step either i or j (or both) decreases
			float x = xRange.mapBackToRange(i, 0, resolution - 1);
			float y = yRange.mapBackToRange(j, 0, resolution - 1);
			prev.add(new XYPoint(x, y));

			if (i <= 0 || j <= 0) {
				return;
			}

			float west = twoDimensionalHistogram[i - 1][j];
			float southwest = twoDimensionalHistogram[i - 1][j - 1];
			float south = twoDimensionalHistogram[i][j - 1];
			float max = Math.max(Math.max(west, south), southwest);

			if (southwest == max || west == south) {
				i--;
				j--;
			} else if (west == max) {
				// We want to go west
				// First check if the last move was south
				// (but only if there was a previous move)
				if (prev.size() >= 2) {
					final XYPoint last = prev.get(prev.size() - 2);
					if (last.x == x) {
						// If the previous points X coord is the same
						// as this point, the last move was to the south.
						// Thus we remove this current point and just move
						// southwest from the last point.
						prev.remove(prev.size() - 1);
					}
				}
				i--;
			} else {
				// We want to go south
				// First check if the last move was west
				// (but only if there was a previous move)
				if (prev.size() >= 2) {
					final XYPoint last = prev.get(prev.size() - 2);
					if (last.y == y) {
						// If the previous points Y coord is the same
						// as this point, the last move was to the west.
						// Thus we remove this current point and just move
						// southwest from the last point.
						prev.remove(prev.size() - 1);
					}
				}
				j--;
			}
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
	
	public static float[][] getStamp(Distribution dist) {
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
					if (stamp[i][j]>0.0f) {
						int localY=yIndex+j-stampRadius;
						if (localY>=0&&localY<twoDimensionalHistogram[i].length) {
							twoDimensionalHistogram[localX][localY]+=stamp[i][j];
						}
					}
				}
			}
		}
	}
}
