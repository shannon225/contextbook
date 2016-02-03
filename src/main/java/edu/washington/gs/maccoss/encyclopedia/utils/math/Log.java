package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class Log {
	private static final double log10=Math.log(10.0);

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
	
	private static final float[] logInts=getLogInts(1000); 
	private static float[] getLogInts(int length) {
		float[] logs=new float[length];
		for (int i=0; i<logs.length; i++) {
			logs[i]=log10(i+1.0f);
		}
		return logs;
	}
	
	/**
	 * an approximation
	 * @param value
	 * @return
	 */
	public static float logFactorial(int value) {
		float sum=0.0f;
		for (int i=0; i<value; i++) {
			sum+=logInts[i];
		}
		return sum;
	}
}
