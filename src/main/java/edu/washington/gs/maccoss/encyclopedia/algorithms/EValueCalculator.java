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
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.procedure.TFloatFloatProcedure;

public class EValueCalculator {
	private static final float FAILURE_STATE = -3f; // (a really bad e-value)

	private final int[] counts=new int[100];

	private final int n;
	private final float m;
	private final float b;

	private final float maxScore;
	private final float maxRT;
	
	private final float minScore;
	private final float binSize;
	

	public EValueCalculator(TFloatFloatHashMap scoreMap, float minScore) {
		this(scoreMap, minScore, estimateBinSize(scoreMap));
	}
	private static float estimateBinSize(TFloatFloatHashMap scoreMap) {

		float[] scoreArray=scoreMap.values();
		float stdev=General.stdev(scoreArray);
		if (scoreArray.length==0||stdev==0.0f) {
			return 1;
		}
		
		// Silverman's (1986) rule of thumb (wikipedia)
		// bandwidth (FWHM for a kernel) for density estimation is similar to binSize in a histogram
		return stdev*(float)Math.pow(4.0/3.0/scoreArray.length, 1.0/5.0);
	}
	/**
	 * 
	 * @param scoreMap keys are RT or index, values are scores
	 * @param minScore
	 * @param binSize
	 */
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
		
		if (scores.size()<3) {
			// use default values if the statistics in the survival function is too meager
			// these are defaults from X!Tandem
			n=scores.size();
			m=-0.25f;
			b=3.5f;
		} else {
			float[] scoreArray = scores.toArray();
			float[] countArray = lnCounts.toArray();
			Pair<Float, Float> equation=LinearRegression.getRegression(scoreArray, countArray);
			n=scores.size();
			
			if (equation.x>=-0.25f) {
				// if the slope is off (or non-negative) then revert back to X!Tandem defaults
				m=-0.25f;
				b=3.5f;
			} else {
				m=equation.x;
				b=equation.y;
			}
		}
	}
	public float getNegLnEValue(float score) {
		if (m>0) return FAILURE_STATE;
		float e=-(score*m+b);
		if (Float.isNaN(e)) {
			return FAILURE_STATE;
		} else if (e<-3) {
			return FAILURE_STATE;
		} else {
			return e;
		}
	}
	
	public float getNegLnEValue() {
		return getNegLnEValue(maxScore);
	}
	public int getN() {
		return n;
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
