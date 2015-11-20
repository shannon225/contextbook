package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class General {
	public static float[] firstDerivative(float[] v) {
		float[] d=new float[v.length-1];
		for (int i=1; i<v.length; i++) {
			d[i-1]=v[i]-v[i-1];
		}
		return d;
	}

	public static float[] protectedSqrt(float[] v) {
		float[] r=new float[v.length];
		for (int i=0; i<r.length; i++) {
			if (v[i]>0) {
				r[i]=(float)Math.sqrt(v[i]);
			}
		}
		return r;
	}
	
	public static float mean(float[] v) {
		float sum=sum(v);
		return sum/v.length;
	}

	public static float sum(float[] v) {
		float sum=0.0f;
		for (int i=0; i<v.length; i++) {
			sum+=v[i];
		}
		return sum;
	}
	
	public static float stdev(float[] v) {
		if (v.length==0) return 0.0f;
		
		float m=mean(v);
		float sumSquares=0.0f;
		for (int i=0; i<v.length; i++) {
			float diff=v[i]-m;
			sumSquares+=diff*diff;
		}
		
		return (float)Math.sqrt(sumSquares/v.length);
	}
	
	public static double mean(double[] v) {
		double sum=sum(v);
		return sum/v.length;
	}

	public static double sum(double[] v) {
		double sum=0.0;
		for (int i=0; i<v.length; i++) {
			sum+=v[i];
		}
		return sum;
	}
	
	public static double stdev(double[] v) {
		if (v.length==0) return 0.0;
		
		double m=mean(v);
		double sumSquares=0.0;
		for (int i=0; i<v.length; i++) {
			double diff=v[i]-m;
			sumSquares+=diff*diff;
		}
		
		return Math.sqrt(sumSquares/v.length);
	}

	public static double[] add(double[] v1, double[] v2) {
		assert(v1.length==v2.length);
		
		double[] r=new double[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]+v2[i];
		}
		return r;
	}
}
