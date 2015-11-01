package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearRegression;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Log;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TFloatFloatHashMap;
import gnu.trove.procedure.TFloatFloatProcedure;

public class EValueCalculator {
	private final int[] counts=new int[100];
	
	private final float m;
	private final float b;
	private final float negLog10EValue;

	float maxScore=-Float.MAX_VALUE;
	float maxRT=0.0f;
	
	public EValueCalculator(TFloatFloatHashMap scoreMap) {
		scoreMap.forEachEntry(new TFloatFloatProcedure() {
			public boolean execute(float arg0, float arg1) {
				if (arg1>maxScore) {
					maxScore=arg1;
					maxRT=arg0;
				}
				return true;
			}
		});
		scoreMap.forEachEntry(new TFloatFloatProcedure() {
			public boolean execute(float arg0, float arg1) {
				int index=getIndex(arg1);
				counts[index]++;
				return true;
			}
		});
		
		TFloatArrayList scores=new TFloatArrayList();
		TFloatArrayList log10Counts=new TFloatArrayList();
		for (int i = 0; i < counts.length; i++) {
			if (counts[i]>1) {
				// don't count singletons, there are a ton and they throw off the extrapolation
				scores.add(getScore(i));
				log10Counts.add(Log.log10(counts[i]));
			}
		}
		
		Pair<Float, Float> equation=LinearRegression.getRegression(scores.toArray(), log10Counts.toArray());
		m=equation.x;
		b=equation.y;
		negLog10EValue=-(maxScore*m+b);
	}
	
	public float getNegLog10EValue() {
		return negLog10EValue;
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
		return (index+0.5f)/counts.length*maxScore;
	}
	
	private int getIndex(float intensity) {
		int index=Math.round(intensity/maxScore*counts.length);
		if (index>=counts.length) {
			index=counts.length-1;
		}
		if (index<0) {
			throw new RuntimeException("No score bins available for scoring system! Empty library spectrum?");
		}
		return index;
	}
}
