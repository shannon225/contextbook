package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;

public class SkylineSGFilter {
	public static float[] paddedSavitzkyGolaySmooth(float[] intRaw) {
		float[] padded=new float[intRaw.length+8];
		System.arraycopy(intRaw, 0, padded, 4, intRaw.length);
		float[] smoothed=savitzkyGolaySmooth(padded);
		float[] unpadded=new float[intRaw.length];
		System.arraycopy(smoothed, 4, unpadded, 0, unpadded.length);
		return unpadded;
	}

	public static float[] savitzkyGolaySmooth(float[] intRaw) {
		if (intRaw==null||intRaw.length<9) return intRaw;
		float[] intSmooth=new float[intRaw.length];
		System.arraycopy(intRaw, 0, intSmooth, 0, 4);
		for (int i=4; i<intRaw.length-4; i++) {
			double sum=59*intRaw[i]+54*(intRaw[i-1]+intRaw[i+1])+39*(intRaw[i-2]+intRaw[i+2])+14*(intRaw[i-3]+intRaw[i+3])-21*(intRaw[i-4]+intRaw[i+4]);
			if (sum<0f) sum=0f;
			intSmooth[i]=(float)(sum/231);
		}
		System.arraycopy(intRaw, intRaw.length-4, intSmooth, intSmooth.length-4, 4);
		return intSmooth;
	}

	public static double[] paddedSavitzkyGolaySmooth(double[] intRaw) {
		double[] padded=new double[intRaw.length+8];
		System.arraycopy(intRaw, 0, padded, 4, intRaw.length);
		double[] smoothed=savitzkyGolaySmooth(padded);
		double[] unpadded=new double[intRaw.length];
		System.arraycopy(smoothed, 4, unpadded, 0, unpadded.length);
		return unpadded;
	}

	public static double[] savitzkyGolaySmooth(double[] intRaw) {
		if (intRaw==null||intRaw.length<9) return intRaw;
		double[] intSmooth=new double[intRaw.length];
		System.arraycopy(intRaw, 0, intSmooth, 0, 4);
		for (int i=4; i<intRaw.length-4; i++) {
			double sum=59*intRaw[i]+54*(intRaw[i-1]+intRaw[i+1])+39*(intRaw[i-2]+intRaw[i+2])+14*(intRaw[i-3]+intRaw[i+3])-21*(intRaw[i-4]+intRaw[i+4]);
			if (sum<0f) sum=0f;
			intSmooth[i]=(sum/231);
		}
		System.arraycopy(intRaw, intRaw.length-4, intSmooth, intSmooth.length-4, 4);
		return intSmooth;
	}

	/**
	 * assumes even incremented X values!
	 * 
	 * @param trace
	 * @return
	 */
	public static XYTrace paddedSavitzkyGolaySmooth(XYTrace trace) {
		Pair<double[], double[]> values=trace.toArrays();
		double[] smoothedY=paddedSavitzkyGolaySmooth(values.y);
		return new XYTrace(values.x, smoothedY, trace.getType(), trace.getName(), trace.getColor().orElse(null), trace.getThickness().orElse(null));
	}
}
