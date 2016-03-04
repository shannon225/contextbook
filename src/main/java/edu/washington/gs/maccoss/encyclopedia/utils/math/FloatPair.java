package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.ArrayList;
import java.util.Comparator;

public class FloatPair implements Comparable<FloatPair> {
	private final float one;
	private final float two;

	public FloatPair(float one, float two) {
		this.one=one;
		this.two=two;
	}

	public float getOne() {
		return one;
	}

	public float getTwo() {
		return two;
	}

	public static class TwosComparator implements Comparator<FloatPair> {
		@Override
		public int compare(FloatPair a, FloatPair b) {
			if (a==null&&b==null) return 0;
			if (a==null) return -1;
			if (b==null) return 1;

			int c=Float.compare(a.getTwo(), b.getTwo());
			if (c!=0) {
				return c;
			}
			return Float.compare(a.getOne(), b.getOne());
		}
	}

	@Override
	public int compareTo(FloatPair o) {
		if (o==null) return 1;
		int c=Float.compare(one, o.getOne());
		if (c!=0) {
			return c;
		}
		return Float.compare(two, o.getTwo());
	}

	public static FloatPair[] toArray(float[] one, float[] two) {
		FloatPair[] pairs=new FloatPair[one.length];
		for (int i=0; i<pairs.length; i++) {
			pairs[i]=new FloatPair(one[i], two[i]);
		}
		return pairs;
	}

	public static float[] getArray(ArrayList<? extends FloatPair> list, boolean getOne) {
		float[] values=new float[list.size()];
		for (int i=0; i<values.length; i++) {
			values[i]=getOne?list.get(i).getOne():list.get(i).getTwo();
		}
		return values;
	}

	public static float[] getArray(FloatPair[] list, boolean getOne) {
		float[] values=new float[list.length];
		for (int i=0; i<values.length; i++) {
			values[i]=getOne?list[i].getOne():list[i].getTwo();
		}
		return values;
	}

	public static double[] getDoubleArray(ArrayList<? extends FloatPair> list, boolean getOne) {
		double[] values=new double[list.size()];
		for (int i=0; i<values.length; i++) {
			values[i]=getOne?list.get(i).getOne():list.get(i).getTwo();
		}
		return values;
	}

	public static double[] getDoubleArray(FloatPair[] list, boolean getOne) {
		double[] values=new double[list.length];
		for (int i=0; i<values.length; i++) {
			values[i]=getOne?list[i].getOne():list[i].getTwo();
		}
		return values;
	}
}
