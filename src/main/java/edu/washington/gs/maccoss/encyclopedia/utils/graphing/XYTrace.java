package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.procedure.TFloatFloatProcedure;

public class XYTrace {
	private final String name;
	private final ArrayList<XYPoint> points;
	private final GraphType type;
	private final Optional<Color> color;
	private final Optional<Float> thickness;
	
	public XYTrace(Spectrum spectrum) {
		color=Optional.empty();
		thickness=Optional.empty();
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

	public XYTrace(Collection<XYPoint> points, GraphType type, String name, Color color, Float thickness) {
		this.color=Optional.ofNullable(color);
		this.thickness=Optional.ofNullable(thickness);
		this.type=type;
		this.points=new ArrayList<XYPoint>(points);
		this.name=name;
		
		Collections.sort(this.points);
	}

	public XYTrace(Collection<XYPoint> points, GraphType type, String name) {
		this(points, type, name, null, null);
	}
	
	public XYTrace(double[] x, double[] y, GraphType type, String name, Color color, Float thickness) {
		this.color=Optional.ofNullable(color);
		this.thickness=Optional.ofNullable(thickness);
		this.type=type;
		this.points=new ArrayList<XYPoint>();
		this.name=name;
		
		assert (x.length==y.length);
		for (int i=0; i<x.length; i++) {
			points.add(new XYPoint(x[i], y[i]));
		}
		Collections.sort(points);
	}
	
	public XYTrace(double[] x, double[] y, GraphType type, String name) {
		this(x, y, type, name, null, null);
	}
	
	public XYTrace(TFloatFloatHashMap map, GraphType type, String name, Color color, Float thickness) {
		this.color=Optional.ofNullable(color);
		this.thickness=Optional.ofNullable(thickness);
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
	public XYTrace(TFloatFloatHashMap map, GraphType type, String name) {
		this(map, type, name, null, null);
	}
	
	public Optional<Color> getColor() {
		return color;
	}
	
	public Optional<Float> getThickness() {
		return thickness;
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

	public static Pair<float[], float[]> toFloatArrays(ArrayList<XYPoint> points) {
		TFloatArrayList xs=new TFloatArrayList();
		TFloatArrayList ys=new TFloatArrayList();
		for (PointInterface point : points) {
			xs.add((float)point.getX());
			ys.add((float)point.getY());
		}
		return new Pair<float[], float[]>(xs.toArray(), ys.toArray());
	}
}
