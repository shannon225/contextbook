package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.procedure.TFloatFloatProcedure;

public class XYTrace {
	private final String name;
	private final ArrayList<XYPoint> points;
	private final GraphType type;
	
	public XYTrace(Spectrum spectrum) {
		this.type=GraphType.spectrum;
		this.points=new ArrayList<XYPoint>();
		this.name=spectrum.getSpectrumName();
		
		double[] mzs=spectrum.getMassArray();
		float[] intensities=spectrum.getIntensityArray();
		
		for (int i=0; i<intensities.length; i++) {
			points.add(new XYPoint(mzs[i], intensities[i]));
		}
		
		Collections.sort(points);
	}

	public XYTrace(Collection<XYPoint> points, GraphType type, String name) {
		this.type=type;
		this.points=new ArrayList<XYPoint>(points);
		this.name=name;
		
		Collections.sort(this.points);
	}
	
	public XYTrace(double[] x, double[] y, GraphType type, String name) {
		this.type=type;
		this.points=new ArrayList<XYPoint>();
		this.name=name;
		
		assert (x.length==y.length);
		for (int i=0; i<x.length; i++) {
			points.add(new XYPoint(x[i], y[i]));
		}
		Collections.sort(points);
	}
	
	public XYTrace(TFloatFloatHashMap map, GraphType type, String name) {
		this.type=type;
		this.points=new ArrayList<XYPoint>();
		this.name=name;

		map.forEachEntry(new TFloatFloatProcedure() {
			public boolean execute(float x, float y) {
				points.add(new XYPoint(x, y));
				return true;
			}
		});
		Collections.sort(points);
	}
	
	public String getName() {
		return name;
	}
	
	public GraphType getType() {
		return type;
	}

	public Pair<double[], double[]> toArrays() {
		return toArrays(points);
	}

	public static Pair<double[], double[]> toArrays(ArrayList<XYPoint> points) {
		TDoubleArrayList xs=new TDoubleArrayList();
		TDoubleArrayList ys=new TDoubleArrayList();
		for (PointInterface point : points) {
			xs.add(point.getX());
			ys.add(point.getY());
		}
		return new Pair<double[], double[]>(xs.toArray(), ys.toArray());
	}
}
