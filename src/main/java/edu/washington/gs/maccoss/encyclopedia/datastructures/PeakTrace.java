package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Optional;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PrecursorIon;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Correlation;
import edu.washington.gs.maccoss.encyclopedia.utils.math.FloatPair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import edu.washington.gs.maccoss.encyclopedia.utils.math.SkylineSGFilter;
import gnu.trove.list.array.TFloatArrayList;

public class PeakTrace <T> implements XYTraceInterface, HasRetentionTime {
	private static final int rtPoints=1000;
	
	private final T ion;
	private final float retentionTimeInSec;
	private final GraphType type;
	private final Color color;
	private final Float thickness;
	private final float[] rt;
	private final float[] intensity;
	private transient PolynomialSplineFunction spline=null;

	public PeakTrace(T ion, float retentionTimeInSec, float[] rt, float[] intensity) {
		this(ion, retentionTimeInSec, GraphType.line, getColor(ion), 1.0f, rt, intensity);
	}
	
	public static Color getColor(Object ion) {
		if (ion instanceof Ion) {
			switch (((PrecursorIon) ion).getIsotope()) {
				case 0: return Color.blue;
				case 1: return Color.magenta;
				case 2: return Color.red;
				case -1: return Color.ORANGE;
				default: return Color.gray;
			}
		}
		return new Color(0, 0, 0, 100);
	}
	
	public PeakTrace(PeakTrace<T> trace, GraphType type, Color color, Float thickness) {
		this((T)trace.ion, trace.retentionTimeInSec, type, color, thickness, trace.rt, trace.intensity);
	}
	
	private PeakTrace(T ion, float retentionTimeInSec, GraphType type, Color color, Float thickness, float[] rt, float[] intensity) {
		this.ion = ion;
		this.retentionTimeInSec=retentionTimeInSec;
		this.type = type;
		this.color = color;
		this.thickness = thickness;
		this.rt = rt;
		this.intensity = intensity;
	}
	
	public PeakTrace<T> updateIntensity(float[] intensity) {
		assert this.rt.length==this.intensity.length:"Number of intensity points does not match number of RT points!";
		
		return new PeakTrace<T>(ion, retentionTimeInSec, type, color, thickness, rt, intensity);
	}
	
	public PeakTrace<T> multiply(float f) {
		return new PeakTrace<T>(ion, retentionTimeInSec, type, color, thickness, rt, General.multiply(intensity, f));
	}
	
	public PeakTrace<T> smooth() {
		float[] smoothed=SkylineSGFilter.paddedSavitzkyGolaySmooth(intensity);
		return new PeakTrace<T>(ion, retentionTimeInSec, type, color, thickness, rt, smoothed);
	}
	
	public T getIon() {
		return ion;
	}
	
	@SuppressWarnings("rawtypes")
	public static PeakTrace<String> getMedianTrace(ArrayList<PeakTrace> traces) {
		if (traces.size()==0) throw new RuntimeException("Trying to generate a median from 0 traces!");
		RetentionTimeBoundary boundary = getBoundary(traces);
		return getMedianTrace(traces, boundary);
	}
	
	@SuppressWarnings("rawtypes")
	public static PeakTrace<String> getMedianTrace(ArrayList<PeakTrace> traces, RetentionTimeBoundary boundary) {
		if (traces.size()==0) throw new RuntimeException("Trying to generate a median from 0 traces!");
		float rtRange=boundary.getMaxRT()-boundary.getMinRT();
		
		TFloatArrayList rtList=new TFloatArrayList();
		TFloatArrayList intensityList=new TFloatArrayList();
		for (int i = 0; i < rtPoints; i++) {
			float thisRT=boundary.getMinRT()+rtRange*(i/(float)rtPoints);
			rtList.add(thisRT);
			TFloatArrayList intensities=new TFloatArrayList();
			for (PeakTrace trace : traces) {
				intensities.add(trace.getIntensity(thisRT));
			}
			intensityList.add(QuickMedian.median(intensities.toArray()));
		}
		
		return new PeakTrace<String>("Median ("+traces.size()+"peaks, ~"+Math.round(boundary.getRetentionTimeInSec()/60)+" mins)", boundary.getRetentionTimeInSec(), GraphType.boldline, Color.red, 3.0f, rtList.toArray(), intensityList.toArray()); 
	}

	@SuppressWarnings("rawtypes")
	private static RetentionTimeBoundary getBoundary(ArrayList<PeakTrace> traces) {
		float min=Float.MAX_VALUE;
		float max=-Float.MAX_VALUE;
		// find the narrowest range
		float averageRT=0.0f;
		for (PeakTrace trace : traces) {
			if (min>trace.getMinRT()) min=trace.getMinRT();
			if (max<trace.getMaxRT()) max=trace.getMaxRT();
			averageRT+=trace.getRetentionTimeInSec();
		}
		averageRT=averageRT/traces.size();
		
		RetentionTimeBoundary boundary=new RetentionTimeBoundary(min, averageRT, max);
		return boundary;
	}
	
	public float getCorrelation(@SuppressWarnings("rawtypes") PeakTrace trace, RetentionTimeBoundary boundary) {
		TFloatArrayList thisIntensities=new TFloatArrayList();
		TFloatArrayList altIntensities=new TFloatArrayList();
		for (int i = 0; i < rt.length; i++) { 
			if (boundary.contains(rt[i])) {
				thisIntensities.add(intensity[i]);
				altIntensities.add(trace.getIntensity(rt[i]));
			}
		}
		return Correlation.getPearsons(thisIntensities.toArray(), altIntensities.toArray());
	}
	
	@Override
	public float getRetentionTimeInSec() {
		return retentionTimeInSec;
	}
	
	public float getIntensity(float retentionTime) {
		if (retentionTime<getMinRT()) return 0.0f;
		if (retentionTime>getMaxRT()) return 0.0f;
		
		if (spline==null) {
			// ok if this gets recalculated, since it's the same every time
			
			if (rt.length==0) return 0; // not enough points to calculate the spline
			
			// force monoisotonic increase in case some scans happen quicker than captured in the mzML floating point value
			float lastRT=-Float.MAX_VALUE;
			for (int i = 0; i < rt.length; i++) {
				if (rt[i]<=lastRT) {
					int bits = Float.floatToIntBits(lastRT)+1;
					rt[i]=Float.intBitsToFloat(bits);
				}
				lastRT=rt[i];
			}
			spline=new SplineInterpolator().interpolate(General.toDoubleArray(rt), General.toDoubleArray(intensity));
		}
		float value = (float)spline.value(retentionTime);
		if (value>0) {
			return value;
		} else {
			return 0;
		}
	}
	
	private static final int numBins=10;
	public float getApexRT(RetentionTimeBoundary boundary) {
		float apex=0.0f;
		int index=-1;
		for (int i = 0; i < rt.length; i++) {
			if (boundary.contains(rt[i])) {
				if (intensity[i]>apex) {
					apex=intensity[i];
					index=i;
				}
			}
		}
		
		int lowerIndex=Math.max(0, index-1);
		int upperIndex=Math.min(rt.length-1, index+1);
		if (index==-1) return 0.0f;
		if (lowerIndex==upperIndex) return intensity[index];
		
		float apexRT = getApex(Math.max(boundary.getMinRT(), rt[lowerIndex]), Math.min(boundary.getMaxRT(), rt[upperIndex]));
		return apexRT;
	}

	// iterate to find target RT until the intensity doesn't change by more than 1%
	private float getApex(float rtStart, float rtStop) {
		float rtRange=rtStop-rtStart;
		float increment=rtRange/numBins;

		float minIntensity=Float.MAX_VALUE;
		float apexIntensity=0.0f;
		float apexRT=0.0f;
		for (int i=0; i<numBins; i++) {
			float localRT=rtStart+i*increment;
			float localIntensity=getIntensity(localRT);
			if (localIntensity>apexIntensity) {
				apexIntensity=localIntensity;
				apexRT=localRT;
			}
			if (localIntensity<minIntensity) {
				minIntensity=localIntensity;
			}
		}
		
		if (minIntensity/apexIntensity<0.99f) {
			return getApex(apexRT-increment, apexRT+increment);
		} else {
			return apexRT;
		}
	}
	
	/**
	 * integrates using actual underlying points, not with spline fitting
	 * @param boundary
	 * @return
	 */
	public float integrateIntensity(RetentionTimeBoundary boundary, boolean backgroundSubtract) {
		FloatPair intensityPair = integrateIntensityAndBackground(boundary);

		if (backgroundSubtract) {
			return Math.max(0.0f, intensityPair.getOne()-intensityPair.getTwo());
		} else {
			return intensityPair.getOne();
		}
	}

	public FloatPair integrateIntensityAndBackground(RetentionTimeBoundary boundary) {
		float totalIntensity=0.0f;
		
		float leftRT=Math.max(rt[0], boundary.getMinRT());
		float leftIntensity=getIntensity(leftRT);
		float rightRT=Math.min(rt[rt.length-1], boundary.getMaxRT());
		float rightIntensity=getIntensity(rightRT);
		float background=calculateTrapezoid(leftRT, rightRT, leftIntensity, rightIntensity);
		
		boolean checkFirstPoint=true;
		int lastValid=-1;
		for (int i = 0; i < rt.length; i++) {
			if (boundary.contains(rt[i])) {
				if (checkFirstPoint) {
					totalIntensity+=calculateTrapezoid(leftRT, rt[i], leftIntensity, intensity[i]);
					checkFirstPoint=false;
				} else {
					totalIntensity+=calculateTrapezoid(rt[i-1], rt[i], intensity[i-1], intensity[i]);
				}
				lastValid=i;
			}
		}
		if (lastValid==-1) return new FloatPair(0.0f, 0.0f);
		
		totalIntensity+=calculateTrapezoid(rt[lastValid], rightRT, intensity[lastValid], rightIntensity);
		
		FloatPair intensityPair=new FloatPair(totalIntensity, background);
		return intensityPair;
	}
	
	protected static float calculateTrapezoid(float rt1, float rt2, float int1, float int2) {
		return (rt2-rt1)*(int1+int2)/2.0f;
	}
	
	@Override
	public Optional<Color> getColor() {
		return Optional.of(color);
	}

	@Override
	public Optional<Float> getThickness() {
		return Optional.of(thickness);
	}

	@Override
	public String getName() {
		return ion.toString();
	}

	@Override
	public GraphType getType() {
		return type;
	}

	@Override
	public Pair<double[], double[]> toArrays() {
		return new Pair<double[], double[]>(General.toDoubleArray(rt), General.toDoubleArray(intensity));
	}

	@Override
	public int size() {
		return rt.length;
	}

	public float getMinRT() {
		return rt[0];
	}
	public float getMaxRT() {
		return rt[rt.length-1];
	}

	public float[] getRt() {
		return rt;
	}

	public float[] getIntensity() {
		return intensity;
	}
}

