package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class Log {
	private static final double log10=Math.log(10.0);
	private static final double log2=Math.log(2.0);

	public static float protectedLn(float v) {
		if (v<=0) return 0.0f;
		return (float)(Math.log(v));
	}

	public static float protectedLog10(float v) {
		if (v<=0) return 0.0f;
		return (float)(Math.log(v)/log10);
	}

	public static double protectedLogLn(double v) {
		if (v<=0) return 0.0;
		return Math.log(v);
	}

	public static double protectedLog10(double v) {
		if (v<=0) return 0.0;
		return Math.log(v)/log10;
	}

	public static double log10(double v) {
		return Math.log(v)/log10;
	}
	public static float log10(float v) {
		return (float)(Math.log(v)/log10);
	}
	public static double log2(double v) {
		return Math.log(v)/log2;
	}
	public static float log2(float v) {
		return (float)(Math.log(v)/log2);
	}

	public static float[] logLn(float[] v) {
		float[] r=new float[v.length];
		for (int i=0; i<r.length; i++) {
			r[i]=(float)Math.log(v[i]);
		}
		return r;
	}

	public static double[] logLn(double[] v) {
		double[] r=new double[v.length];
		for (int i=0; i<r.length; i++) {
			r[i]=Math.log(v[i]);
		}
		return r;
	}

	public static float[] log10(float[] v) {
		float[] r=new float[v.length];
		for (int i=0; i<r.length; i++) {
			r[i]=log10(v[i]);
		}
		return r;
	}

	public static double[] log10(double[] v) {
		double[] r=new double[v.length];
		for (int i=0; i<r.length; i++) {
			r[i]=log10(v[i]);
		}
		return r;
	}
	
	private static final float[] logFacts=getLogFacts(1000); 
	private static float[] getLogFacts(int length) {
		float[] logs=new float[length];
		float sum=0.0f;
		for (int i=0; i<logs.length-1; i++) {
			sum+=log10(i+1.0f);
			logs[i+1]=sum;
		}
		return logs;
	}
	
	/**
	 * an approximation
	 * @param value
	 * @return
	 */
	public static float logFactorial(int value) {
		return logFacts[value];
	}
}
