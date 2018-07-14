package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.List;
import java.util.Set;

import edu.washington.gs.maccoss.encyclopedia.datastructures.IntRange;

public class General {
	public static String formatCellToWidth(String s, int w) {
		char[] ca=new char[w];
		for (int i=0; i<ca.length; i++) {
			ca[i]=' ';
		}
		for (int i=0; i<s.length(); i++) {
			if (i<ca.length-1) {
				ca[i]=s.charAt(i);
			}
		}
		return new String(ca);
	}
	
	public static int numberOfOccurances(String s, String target) {
		int indexOf=s.indexOf(target);
		int count=0;
		while (indexOf>=0) {
			count++;
			indexOf=s.indexOf(target, indexOf+target.length());
		}
		return count;
	}

	public static double[][] transposeMatrix(double[][] m) {
		double[][] temp = new double[m[0].length][m.length];
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[0].length; j++)
				temp[j][i] = m[i][j];
		return temp;
	}
	
	public static float[] toFloatArray(Float[] a) {
		float[] f=new float[a.length];
		for (int i=0; i<f.length; i++) {
			f[i]=a[i]==null?0.0f:a[i].floatValue();
		}
		return f;
	}
	
	public static float[] toFloatArray(double[] a) {
		float[] f=new float[a.length];
		for (int i=0; i<f.length; i++) {
			f[i]=(float)a[i];
		}
		return f;
	}
	
	public static double[] toDoubleArray(float[] a) {
		double[] f=new double[a.length];
		for (int i=0; i<f.length; i++) {
			f[i]=a[i];
		}
		return f;
	}
	
	public static boolean equals(int[] a, int[] b) {
		if (a.length!=b.length) return false;
		for (int i=0; i<b.length; i++) {
			if (a[i]!=b[i]) return false;
		}
		return true;
	}
	
	public static String[] concatenate(String[] a, String... s) {
		String[] r=new String[a.length+s.length];
		System.arraycopy(a, 0, r, 0, a.length);
		System.arraycopy(s, 0, r, a.length, s.length);
		return r;
	}

	public static String[] concatenate(String[]... a) {
		int length=0;
		for (int i=0; i<a.length; i++) {
			length+=a[i].length;
		}
		String[] r=new String[length];
		
		int lastIndex=0;
		for (int i=0; i<a.length; i++) {
			System.arraycopy(a[i], 0, r, lastIndex, a[i].length);
			lastIndex+=a[i].length;
		}
		return r;
	}
	
	public static float[] concatenate(float[] a, float... f) {
		float[] r=new float[a.length+f.length];
		System.arraycopy(a, 0, r, 0, a.length);
		System.arraycopy(f, 0, r, a.length, f.length);
		return r;
	}

	public static float[] concatenate(float[]... a) {
		int length=0;
		for (int i=0; i<a.length; i++) {
			length+=a[i].length;
		}
		float[] r=new float[length];
		
		int lastIndex=0;
		for (int i=0; i<a.length; i++) {
			System.arraycopy(a[i], 0, r, lastIndex, a[i].length);
			lastIndex+=a[i].length;
		}
		return r;
	}
	
	public static String toString(double[] i, String delim) {
		StringBuilder sb=new StringBuilder();
		for (Object g : i) {
			if (sb.length()>0) {
				sb.append(delim);
			}
			sb.append(g);
		}
		return sb.toString();
	}
	
	public static String toString(Object[] i, String delim) {
		StringBuilder sb=new StringBuilder();
		for (Object g : i) {
			if (sb.length()>0) {
				sb.append(delim);
			}
			sb.append(g);
		}
		return sb.toString();
	}
	
	public static String toString(Object[] i) {
		return toString(i, ",");
	}
	
	public static String toString(@SuppressWarnings("rawtypes") List i) {
		StringBuilder sb=new StringBuilder();
		for (Object g : i) {
			if (sb.length()>0) {
				sb.append(',');
			}
			sb.append(g);
		}
		return sb.toString();
	}
	
	public static String toString(@SuppressWarnings("rawtypes") Set i) {
		if (i==null) return null;
		
		StringBuilder sb=new StringBuilder();
		for (Object g : i) {
			if (sb.length()>0) {
				sb.append(',');
			}
			sb.append(g);
		}
		return sb.toString();
	}
	
	public static String toString(int[] i) {
		if (i==null) return null;
		
		StringBuilder sb=new StringBuilder();
		for (int g : i) {
			if (sb.length()>0) {
				sb.append(',');
			}
			sb.append(g);
		}
		return sb.toString();
	}
	
	public static String toString(byte[] i) {
		if (i==null) return null;
		
		StringBuilder sb=new StringBuilder();
		for (int g : i) {
			if (sb.length()>0) {
				sb.append(',');
			}
			sb.append(g);
		}
		return sb.toString();
	}
	
	public static String toString(float[] f, String delim) {
		if (f==null) return null;
		
		StringBuilder sb=new StringBuilder();
		for (float g : f) {
			if (sb.length()>0) {
				sb.append(delim);
			}
			sb.append(g);
		}
		return sb.toString();
	}
	
	public static String toString(float[] f) {
		return toString(f, ",");
	}
	
	public static String toString(double[] i) {
		if (i==null) return null;
		
		StringBuilder sb=new StringBuilder();
		for (Object g : i) {
			if (sb.length()>0) {
				sb.append(',');
			}
			sb.append(g);
		}
		return sb.toString();
	}
	
	public static float[] normalizeAndBackgroundSubtract(float[] v, IntRange range) {
		v=v.clone();
		
		int stop=Math.min(v.length-1, range.getStop()+1);
		int start=Math.max(0, range.getStart()-1);
		float deltaY=v[stop]-v[start];
		float deltaX=stop-start;
		if (deltaX==0.0f) return normalize(v, range); //new float[v.length];
		
		float max=General.max(extract(v, range));
		if (v[start]>=max||v[stop]>=max) return normalize(v, range);
		
		float m=deltaY/deltaX;
		float b=v[stop]-m*stop;
		
		for (int i=0; i<v.length; i++) {
			if (range.contains(i)) {
				float background=m*i+b;
				if (background>v[i]) {
					v[i]=0.0f;
				} else if (background>0.0f) {
					v[i]=v[i]-background;
				}
			}
		}
		
		return normalize(v, range);
	}
	
	public static float[] normalize(float[] v, IntRange range) {
		float sum=sum(v, range);
		if (sum==0.0f) {
			return new float[v.length];
		}
		return divide(v, sum);
	}
	
	public static int[] extract(int[] v, IntRange range) {
		int[] r=new int[range.getRange()];
		for (int i=0; i<r.length; i++) {
			r[i]=v[i+range.getStart()];
		}
		return r;
	}
	
	public static float[] extract(float[] v, IntRange range) {
		float[] r=new float[range.getRange()];
		for (int i=0; i<r.length; i++) {
			r[i]=v[i+range.getStart()];
		}
		return r;
	}
	
	public static float[] normalize(float[] v) {
		return normalize(v, new IntRange(0, v.length-1));
	}
	
	public static double[] divide(double[] v, double d) {
		double[] f=new double[v.length];
		for (int i=0; i<v.length; i++) {
			f[i]=v[i]/d;
		}
		return f;
	}
	
	public static float[] divide(float[] v, float d) {
		float[] f=new float[v.length];
		for (int i=0; i<v.length; i++) {
			f[i]=v[i]/d;
		}
		return f;
	}
	
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
	
	public static float mean(int[] v) {
		if (v.length==0) return 0.0f;
		
		float sum=sum(v);
		return sum/(float)v.length;
	}
	
	public static float mean(float[] v) {
		if (v.length==0) return 0.0f;
		
		float sum=sum(v);
		return sum/v.length;
	}
	
	public static float mean(float[][] v) {
		if (v.length==0) return 0.0f;
		int length=0;
		float sum=0.0f;
		for (int i=0; i<v.length; i++) {
			length+=v[i].length;
			sum+=sum(v[i]);
		}
		return sum/length;
	}
	
	public static float mean(float[] v, int startIndex, int stopIndex) {
		float sum=0.0f;
		int count=0;
		for (int i=startIndex; i<=stopIndex; i++) {
			sum+=v[i];
			count++;
		}
		if (count==0) return 0.0f;
		return sum/count;
	}

	public static float sum(float[] v, IntRange range) {
		if (v==null||v.length==0) return 0.0f;
		
		float sum=0.0f;
		for (int i=range.getStart(); i<=range.getStop(); i++) {
			sum+=v[i];
		}
		return sum;
	}

	public static float sum(float[] v) {
		if (v==null||v.length==0) return 0.0f;
		return sum(v, new IntRange(0, v.length-1));
	}

	public static int sum(int[] v) {
		int sum=0;
		for (int i=0; i<v.length; i++) {
			sum+=v[i];
		}
		return sum;
	}
	
	public static float variance(float[] v) {
		if (v.length==0) return 0.0f;
		
		float m=mean(v);
		float sumSquares=0.0f;
		for (int i=0; i<v.length; i++) {
			float diff=v[i]-m;
			sumSquares+=diff*diff;
		}
		
		return sumSquares/v.length;
	}
	
	public static float stdev(float[] v) {
		if (v.length==0) return 0.0f;
		
		return (float)Math.sqrt(variance(v));
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

	public static float[] add(float[] v1, float[] v2) {
		assert(v1.length==v2.length);
		
		float[] r=new float[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]+v2[i];
		}
		return r;
	}

	public static double[] multiply(double[] v1, double m) {
		double[] r=new double[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]*m;
		}
		return r;
	}

	public static float[] multiply(float[] v1, float m) {
		float[] r=new float[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]*m;
		}
		return r;
	}

	public static float[] multiply(float[] v1, float[] v2) {
		assert(v1.length==v2.length);
		
		float[] r=new float[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]*v2[i];
		}
		return r;
	}

	public static float[] add(float[] v1, float v) {
		float[] r=new float[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]+v;
		}
		return r;
	}

	public static float[] subtract(float[] v1, float v) {
		float[] r=new float[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]-v;
		}
		return r;
	}

	public static double[] subtract(double[] v1, double v) {
		double[] r=new double[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]-v;
		}
		return r;
	}

	public static double[] subtract(double[] v1, double[] v2) {
		double[] r=new double[v1.length];
		for (int i=0; i<r.length; i++) {
			r[i]=v1[i]-v2[i];
		}
		return r;
	}
	
	public static double max(double[] v) {
		double max=-Double.MAX_VALUE;
		for (int i=0; i<v.length; i++) {
			if (v[i]>max) {
				max=v[i];
			}
		}
		return max;
	}
	
	public static float max(float[] v) {
		float max=-Float.MAX_VALUE;
		for (int i=0; i<v.length; i++) {
			if (v[i]>max) {
				max=v[i];
			}
		}
		return max;
	}
	
	public static int max(int[] v) {
		int max=-Integer.MAX_VALUE;
		for (int i=0; i<v.length; i++) {
			if (v[i]>max) {
				max=v[i];
			}
		}
		return max;
	}
	
	public static double min(double[] v) {
		double min=Double.MAX_VALUE;
		for (int i=0; i<v.length; i++) {
			if (v[i]<min) {
				min=v[i];
			}
		}
		return min;
	}
	
	public static float min(float[] v) {
		float min=Float.MAX_VALUE;
		for (int i=0; i<v.length; i++) {
			if (v[i]<min) {
				min=v[i];
			}
		}
		return min;
	}
	
	public static int min(int[] v) {
		int min=Integer.MAX_VALUE;
		for (int i=0; i<v.length; i++) {
			if (v[i]<min) {
				min=v[i];
			}
		}
		return min;
	}
}
