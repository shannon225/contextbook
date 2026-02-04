package edu.washington.gs.maccoss.encyclopedia.utils.math;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;

public class SkylineSGFilter {
	public static float[] paddedSavitzkyGolaySmooth(float[] intRaw) {
		return paddedSavitzkyGolay9PointSmooth(intRaw);
	}
	public static float[] savitzkyGolaySmooth(float[] rawArray) {
		return savitzkyGolay9PointSmooth(rawArray);
	}

	public static XYTrace paddedSavitzkyGolaySmooth(XYTraceInterface trace) {
		return paddedSavitzkyGolay9PointSmooth(trace);
	}
	
	public static float[] paddedSavitzkyGolay9PointSmooth(float[] intRaw) {
		float[] padded=new float[intRaw.length+8];
		System.arraycopy(intRaw, 0, padded, 4, intRaw.length);
		float[] smoothed=savitzkyGolay9PointSmooth(padded);
		float[] unpadded=new float[intRaw.length];
		System.arraycopy(smoothed, 4, unpadded, 0, unpadded.length);
		return unpadded;
	}

	public static float[] savitzkyGolay9PointSmooth(float[] rawArray) {
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
	public static float[] paddedSavitzkyGolay7PointSmooth(float[] intRaw) {
		float[] padded=new float[intRaw.length+6];
		System.arraycopy(intRaw, 0, padded, 3, intRaw.length);
		float[] smoothed=savitzkyGolay7PointSmooth(padded);
		float[] unpadded=new float[intRaw.length];
		System.arraycopy(smoothed, 3, unpadded, 0, unpadded.length);
		return unpadded;
	}

	public static float[] savitzkyGolay7PointSmooth(float[] rawArray) {
		if (rawArray==null||rawArray.length<7) return rawArray;
		float[] smoothedArray=new float[rawArray.length];
		for (int i = 0; i < 3; i++) {
			smoothedArray[i]=rawArray[i];
		}
		for (int i = rawArray.length-3; i < rawArray.length; i++) {
			smoothedArray[i]=rawArray[i];
		}
		for (int i=3; i<rawArray.length-3; i++) {
			float sum=7f*rawArray[i]
					+6f*(rawArray[i-1]+rawArray[i+1])
					+3f*(rawArray[i-2]+rawArray[i+2])
					-2f*(rawArray[i-3]+rawArray[i+3]);
			
			if (sum<0f) sum=0f;
			smoothedArray[i]=sum/21;
		}
		return smoothedArray;
	}
	public static float[] paddedSavitzkyGolay5PointSmooth(float[] intRaw) {
		float[] padded=new float[intRaw.length+4];
		System.arraycopy(intRaw, 0, padded, 2, intRaw.length);
		float[] smoothed=savitzkyGolay5PointSmooth(padded);
		float[] unpadded=new float[intRaw.length];
		System.arraycopy(smoothed, 2, unpadded, 0, unpadded.length);
		return unpadded;
	}

	public static float[] savitzkyGolay5PointSmooth(float[] rawArray) {
		if (rawArray==null||rawArray.length<5) return rawArray;
		float[] smoothedArray=new float[rawArray.length];
		for (int i = 0; i < 2; i++) {
			smoothedArray[i]=rawArray[i];
		}
		for (int i = rawArray.length-2; i < rawArray.length; i++) {
			smoothedArray[i]=rawArray[i];
		}
		for (int i=2; i<rawArray.length-2; i++) {
			float sum=17f*rawArray[i]
					+12f*(rawArray[i-1]+rawArray[i+1])
					-3f*(rawArray[i-2]+rawArray[i+2]);
			
			if (sum<0f) sum=0f;
			smoothedArray[i]=sum/35;
		}
		return smoothedArray;
	}
	public static float[] paddedSavitzkyGolay3PointSmooth(float[] intRaw) {
		float[] padded=new float[intRaw.length+2];
		System.arraycopy(intRaw, 0, padded, 1, intRaw.length);
		float[] smoothed=savitzkyGolay3PointSmooth(padded);
		float[] unpadded=new float[intRaw.length];
		System.arraycopy(smoothed, 1, unpadded, 0, unpadded.length);
		return unpadded;
	}

	public static float[] savitzkyGolay3PointSmooth(float[] rawArray) {
		if (rawArray==null||rawArray.length<3) return rawArray;
		float[] smoothedArray=new float[rawArray.length];
		for (int i = 0; i < 1; i++) {
			smoothedArray[i]=rawArray[i];
		}
		for (int i = rawArray.length-1; i < rawArray.length; i++) {
			smoothedArray[i]=rawArray[i];
		}
		for (int i=1; i<rawArray.length-1; i++) {
			float sum=rawArray[i]+rawArray[i-1]+rawArray[i+1];
			
			if (sum<0f) sum=0f;
			smoothedArray[i]=sum/3;
		}
		return smoothedArray;
	}

	/**
	 * assumes even incremented X values!
	 * 
	 * @param trace
	 * @return
	 */
	public static XYTrace paddedSavitzkyGolay9PointSmooth(XYTraceInterface trace) {
		Pair<double[], double[]> values=trace.toArrays();
		double[] smoothedY=General.toDoubleArray(paddedSavitzkyGolay9PointSmooth(General.toFloatArray(values.y)));
		return new XYTrace(values.x, smoothedY, trace.getType(), trace.getName(), trace.getColor().orElse(null), trace.getThickness().orElse(null));
	}
	public static XYTrace paddedSavitzkyGolay7PointSmooth(XYTraceInterface trace) {
		Pair<double[], double[]> values=trace.toArrays();
		double[] smoothedY=General.toDoubleArray(paddedSavitzkyGolay7PointSmooth(General.toFloatArray(values.y)));
		return new XYTrace(values.x, smoothedY, trace.getType(), trace.getName(), trace.getColor().orElse(null), trace.getThickness().orElse(null));
	}
	public static XYTrace paddedSavitzkyGolay5PointSmooth(XYTraceInterface trace) {
		Pair<double[], double[]> values=trace.toArrays();
		double[] smoothedY=General.toDoubleArray(paddedSavitzkyGolay5PointSmooth(General.toFloatArray(values.y)));
		return new XYTrace(values.x, smoothedY, trace.getType(), trace.getName(), trace.getColor().orElse(null), trace.getThickness().orElse(null));
	}
	public static XYTrace paddedSavitzkyGolay3PointSmooth(XYTraceInterface trace) {
		Pair<double[], double[]> values=trace.toArrays();
		double[] smoothedY=General.toDoubleArray(paddedSavitzkyGolay3PointSmooth(General.toFloatArray(values.y)));
		return new XYTrace(values.x, smoothedY, trace.getType(), trace.getName(), trace.getColor().orElse(null), trace.getThickness().orElse(null));
	}

	public static XYTrace adjustableSavitzkyGolaySmooth(XYTraceInterface trace) {
		Pair<double[], double[]> values=trace.toArrays();
		float[] yArray=General.toFloatArray(values.y);
		
		float minValue=0.05f*General.max(yArray);
		int count=0;
		for (int i=0; i<yArray.length; i++) {
			if (yArray[i]>minValue) count++;
		}
		
		double[] smoothedY;
		if (count<=1) {
			smoothedY=values.y; 
		} else if (count<3) {
			smoothedY=General.toDoubleArray(paddedSavitzkyGolay3PointSmooth(yArray));
		} else if (count<5) {
			smoothedY=General.toDoubleArray(paddedSavitzkyGolay5PointSmooth(yArray));
		} else if (count<7) {
			smoothedY=General.toDoubleArray(paddedSavitzkyGolay7PointSmooth(yArray));
		} else {
			smoothedY=General.toDoubleArray(paddedSavitzkyGolay9PointSmooth(yArray));
		}
		
		return new XYTrace(values.x, smoothedY, trace.getType(), trace.getName(), trace.getColor().orElse(null), trace.getThickness().orElse(null));
	}
}
