package edu.washington.gs.maccoss.encyclopedia.utils.graphing;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.procedure.TFloatFloatProcedure;

public class XYTrace implements XYTraceInterface {
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
	
	public static double getMaxY(XYTraceInterface[] traces) {
		double max=-Double.MAX_VALUE;
		for (XYTraceInterface xyTrace : traces) {
			for (double y : xyTrace.toArrays().y) {
				if (y>max) {
					max=y;
				}
			}
		}
		return max;
	}
	
	public static double getMaxY(Collection<XYTrace> traces) {
		double max=-Double.MAX_VALUE;
		for (XYTrace xyTrace : traces) {
			double newMax=xyTrace.getMaxY();
			if (newMax>max) {
				max=newMax;
			}
		}
		return max;
	}
	
	public double getMaxYInRange(Range xrange) {
		double max=-Double.MAX_VALUE;
		for (XYPoint xy : points) {
			if (xrange.contains(xy.getX())) {
				if (xy.y>max) {
					max=xy.y;
				}
			}
		}
		return max;
	}
	
	public double getMaxY() {
		double max=-Double.MAX_VALUE;
		for (XYPoint xy : points) {
			if (xy.y>max) {
				max=xy.y;
			}
		}
		return max;
	}
	
	public XYTraceInterface rescaleX(float rescaleX) {
		Pair<double[], double[]> trace=toArrays(points);
		double[] newx=General.multiply(trace.x, rescaleX);
		return new XYTrace(newx, trace.y, type, name, color.orElse(null), thickness.orElse(null));
	}
	
	public XYTrace(double[] x, double[] y, GraphType type, String name, Optional<Color> color, Optional<Float> thickness) {
		this(x, y, type, name, color.orElse(null), thickness.orElse(null));
	}
	
	public XYTrace(ArrayList<XYPoint> points, GraphType type, String name, Optional<Color> color, Optional<Float> thickness) {
		this(points, type, name, color.orElse(null), thickness.orElse(null));
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
		this(x, y, type, name, Optional.ofNullable((Color)null), Optional.ofNullable((Float)null));
	}
	
	public XYTrace(float[] x, float[] y, GraphType type, String name) {
		this(General.toDoubleArray(x), General.toDoubleArray(y), type, name, Optional.ofNullable((Color)null), Optional.ofNullable((Float)null));
	}
	
	public XYTrace(float[] x, float[] y, GraphType type, String name, Color color, Float thickness) {
		this(General.toDoubleArray(x), General.toDoubleArray(y), type, name, color, thickness);
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
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface#getColor()
	 */
	@Override
	public Optional<Color> getColor() {
		return color;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface#getThickness()
	 */
	@Override
	public Optional<Float> getThickness() {
		return thickness;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface#getName()
	 */
	@Override
	public String getName() {
		return name;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface#getType()
	 */
	@Override
	public GraphType getType() {
		return type;
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface#toArrays()
	 */
	@Override
	public Pair<double[], double[]> toArrays() {
		return toArrays(points);
	}
	
	public String toString() {
		Pair<double[], double[]> pair=toArrays(points);
		StringBuilder sb=new StringBuilder("// "+getName()+"\n");
		sb.append("float[] x=new float[] {");
		boolean first=true;
		for (double d : pair.x) {
			if (first) {
				first=false;
			} else {
				sb.append(',');
			}
			sb.append(d);
			sb.append('f');
		}
		sb.append("};\n");
		sb.append("float[] y=new float[] {");
		first=true;
		for (double d : pair.y) {
			if (first) {
				first=false;
			} else {
				sb.append(',');
			}
			sb.append(d);
			sb.append('f');
		}
		sb.append("};\n");
		return sb.toString();
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
