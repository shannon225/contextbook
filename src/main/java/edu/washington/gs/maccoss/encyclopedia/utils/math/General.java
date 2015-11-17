package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class General {

	public static float[] protectedSqrt(float[] v) {
		float[] r=new float[v.length];
		for (int i=0; i<r.length; i++) {
			if (v[i]>0) {
				r[i]=(float)Math.sqrt(v[i]);
			}
		}
		return r;
	}
	
	public static double mean(double[] v) {
		double sum=0.0;
		for (int i=0; i<v.length; i++) {
			sum+=v[i];
		}
		return sum/v.length;
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
}
