package edu.washington.gs.maccoss.encyclopedia.utils.math;

public class Log {
	private static final double log10=Math.log(10.0);

	public static double log10(double v) {
		return Math.log(v)/log10;
	}
	public static float log10(float v) {
		return (float)(Math.log(v)/log10);
	}
}
