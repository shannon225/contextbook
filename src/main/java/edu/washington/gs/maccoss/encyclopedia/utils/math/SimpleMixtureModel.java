package edu.washington.gs.maccoss.encyclopedia.utils.math;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.math.distributions.Distribution;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class SimpleMixtureModel implements RTProbabilityModel {
	private static final int resolution=100; 
	private Distribution positive;
	private Distribution negative;
	private Range range;
	private PolynomialSplineFunction spline;
	
	public SimpleMixtureModel(Distribution positive, Distribution negative, Range range) {
		this.positive=positive;
		this.negative=negative;
		
		TDoubleArrayList xs=new TDoubleArrayList();
		TFloatArrayList probs=new TFloatArrayList();
		
		for (int i=0; i<resolution; i++) {
			float x=(i/(float)resolution)*range.getRange()+range.getStart();
			double p=positive.getProbability(x);
			double n=negative.getProbability(x);
			double sum=p+n;
			double y;
			if (sum==0) {
				y=0.0;
			} else {
				y=p/sum;
			}
			
			xs.add(x);
			probs.add((float)y);
		}
		this.range=new Range(xs.get(0), xs.get(xs.size()-1));
		
		float[] array = probs.toArray();
		array=BackgroundSubtractionFilter.fastMovingMedian(array, resolution/10);
		
		spline=new SplineInterpolator().interpolate(xs.toArray(), General.toDoubleArray(array));
	}

	@Override
	public float getProbability(float retentionTime, float delta) {
		if (range.contains(delta)) {
			return (float)spline.value(delta);
		}
		return 0.0f;
	}
	
	public Distribution getPositive() {
		return positive;
	}
	
	public Distribution getNegative() {
		return negative;
	}
}
