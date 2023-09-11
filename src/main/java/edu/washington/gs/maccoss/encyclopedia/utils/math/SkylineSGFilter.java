package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;

public class SkylineSGFilter {
	public static float[] paddedSavitzkyGolaySmooth(float[] intRaw) {
		float[] padded=new float[intRaw.length+8];
		System.arraycopy(intRaw, 0, padded, 4, intRaw.length);
		float[] smoothed=savitzkyGolaySmooth(padded);
		float[] unpadded=new float[intRaw.length];
		System.arraycopy(smoothed, 4, unpadded, 0, unpadded.length);
		return unpadded;
	}

	public static float[] savitzkyGolaySmooth(float[] rawArray) {
		if (rawArray==null||rawArray.length<9) return rawArray;
		float[] smoothedArray=new float[rawArray.length];
		for (int i = 0; i < 4; i++) {
			smoothedArray[i]=rawArray[i];
		}
		for (int i = rawArray.length-4; i < rawArray.length; i++) {
			smoothedArray[i]=rawArray[i];
		}
		for (int i=4; i<rawArray.length-4; i++) {
			float sum=59f*rawArray[i]
					+54f*(rawArray[i-1]+rawArray[i+1])
					+39f*(rawArray[i-2]+rawArray[i+2])
					+14f*(rawArray[i-3]+rawArray[i+3])
					-21f*(rawArray[i-4]+rawArray[i+4]);
			
			if (sum<0f) sum=0f;
			smoothedArray[i]=sum/231;
		}
		return smoothedArray;
	}

	/**
	 * assumes even incremented X values!
	 * 
	 * @param trace
	 * @return
	 */
	public static XYTrace paddedSavitzkyGolaySmooth(XYTraceInterface trace) {
		Pair<double[], double[]> values=trace.toArrays();
		double[] smoothedY=General.toDoubleArray(paddedSavitzkyGolaySmooth(General.toFloatArray(values.y)));
		return new XYTrace(values.x, smoothedY, trace.getType(), trace.getName(), trace.getColor().orElse(null), trace.getThickness().orElse(null));
	}
}
