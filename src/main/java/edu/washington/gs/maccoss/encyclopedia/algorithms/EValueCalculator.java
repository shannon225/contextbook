package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearRegression;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.procedure.TFloatFloatProcedure;

public class EValueCalculator {
	private final int[] counts=new int[100];
	
	private final float m;
	private final float b;
	private final float neglnEValue;

	private final float maxScore;
	private final float maxRT;
	
	private final float minScore;
	private final float binSize;
	
	public EValueCalculator(TFloatFloatHashMap scoreMap, float minScore, float binSize) {
		this.minScore=minScore;
		this.binSize=binSize;
		
		final float[] intermediateMaxScores=new float[2];
		intermediateMaxScores[0]=-Float.MAX_VALUE;
		intermediateMaxScores[1]=0.0f;
		scoreMap.forEachEntry(new TFloatFloatProcedure() {
			public boolean execute(float arg0, float arg1) {
				if (arg1>intermediateMaxScores[0]) {
					intermediateMaxScores[0]=arg1;
					intermediateMaxScores[1]=arg0;
				} else if (arg1==intermediateMaxScores[0]&&arg0>intermediateMaxScores[1]) {
					// force a repeatable sort order for exact score matches
					intermediateMaxScores[0]=arg1;
					intermediateMaxScores[1]=arg0;
				}
				return true;
			}
		});
		maxScore=intermediateMaxScores[0];
		maxRT=intermediateMaxScores[1];
		
		scoreMap.forEachEntry(new TFloatFloatProcedure() {
			public boolean execute(float arg0, float arg1) {
				if (arg1>=minScore) { // require scores of at least minScore
					int index=Math.round((arg1-minScore)/binSize);
					if (index>=counts.length) {
						index=counts.length-1;
					}
					if (index>=0) {
						counts[index]++;
					}
				}
				return true;
			}
		});
		
		int totalCounts=General.sum(counts);
		int target=Math.round(totalCounts/2f);
		
		int targetIndex=0;
		int currentTotal=0;
		for (int i = counts.length-1; i>=0; i--) {
			currentTotal+=counts[i];
			if (currentTotal>target) {
				targetIndex=i;
				break;
			}
		}
		
		TFloatArrayList scores=new TFloatArrayList();
		TFloatArrayList lnCounts=new TFloatArrayList();
		for (int i = targetIndex; i < counts.length; i++) {
			if (counts[i]>3) {
				// don't count singletons, there are a ton and they throw off the extrapolation
				scores.add(getScore(i));
				lnCounts.add((float)Math.log(counts[i]));
				
				//System.out.println(getScore(i)+"\t"+Math.log(counts[i]));
			}
		}
		
		Pair<Float, Float> equation=LinearRegression.getRegression(scores.toArray(), lnCounts.toArray());
		m=equation.x;
		b=equation.y;
		if (m>0) {
			neglnEValue=-3;
		} else {
			float e = -(maxScore*m+b);
			neglnEValue=e<-3?-3:e;
		}
	}
	public float getNegLnEValue(float score) {
		return -(score*m+b);
	}
	
	public float getNegLnEValue() {
		return neglnEValue;
	}
	public float getB() {
		return b;
	}
	public float getM() {
		return m;
	}

	public float getMaxRT() {
		return maxRT;
	}

	public float getMaxRawScore() {
		return maxScore;
	}
	
	/**
	 * returns the score for the middle of the bin
	 * @param index
	 * @return
	 */
	private float getScore(int index) {
		return (index+0.5f)*binSize+minScore;
	}
	
	public XYTraceInterface[] toTraces() {
		ArrayList<XYPoint> points1=new ArrayList<XYPoint>();
		ArrayList<XYPoint> points2=new ArrayList<XYPoint>();
		for (int i=0; i<counts.length; i++) {
			float intensity=getScore(i);
			points1.add(new XYPoint(intensity, Log.protectedLn(counts[i])));
			points2.add(new XYPoint(intensity, (intensity*m+b)));
		}
		return new XYTraceInterface[] {new XYTrace(points2, GraphType.line, "fit"), new XYTrace(points1, GraphType.area, "histogram")};
	}
	
	public XYTraceInterface[] toUnloggedTraces() {
		ArrayList<XYPoint> points1=new ArrayList<XYPoint>();
		ArrayList<XYPoint> points2=new ArrayList<XYPoint>();
		for (int i=0; i<counts.length; i++) {
			float intensity=getScore(i);
			points1.add(new XYPoint(intensity, counts[i]));
			points2.add(new XYPoint(intensity, Math.pow(Math.E, (intensity*m+b))));
		}
		return new XYTraceInterface[] {new XYTrace(points2, GraphType.line, "fit"), new XYTrace(points1, GraphType.area, "histogram")};
	}
}
