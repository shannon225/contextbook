package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.awt.Color;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class RandomGenerator {
	private static final char[] aas="ACDEFGHIKLMNPQRSTVWY".toCharArray(); 
	/**
	 * rand() from ANSI C, should be faster (and easier to control) than java.util.Random, which uses Longs
	 * @param seed (Please don't use 0 as a seed!)
	 * @return
	 */
	public static int randomInt(int seed) {
		seed=seed*1103515245+12345;
		return seed%2147483647;
	}
	public static float random(int seed) {
		return floatFromRandomInt(randomInt(seed));
	}
	public static float floatFromRandomInt(int random) {
		return ((random/(float)2147483647)+1f)/2f;
	}
	public static int randomIndex(int length, int seed) {
		return (int)(length*floatFromRandomInt(randomInt(seed)));
	}

	public static void shuffle(List<?> list, int seed) {
		Collections.shuffle(list, new Random(seed));
    }
	
	public static String randomSequence(int seed) {
		int length=(int)(20*random(seed))+6;
		
		StringBuilder sb=new StringBuilder();
		for (int i=0; i<length; i++) {
			seed=randomInt(seed);
			char aa=aas[randomIndex(aas.length, seed)];
			sb.append(aa);
		}
		return sb.toString();
	}
	
	public static Color randomColor(int seed) {
		int random=randomInt(seed);
		random=randomInt(random);
		random=randomInt(random);
		random=randomInt(random);
		random=randomInt(random);
		Color c=new Color(random);
		
		float intensifier=(Math.abs(random)%60)/100f;
		float detensifier=(Math.abs(random)%40)/100f;
		
		float red=c.getRed()/255f;
		float green=c.getGreen()/255f;
		float blue=c.getBlue()/255f;
		float[] values=new float[] {red, green, blue};
		float median=QuickMedian.median(values);
		
		if (red>median) {
			red=red+intensifier;
		} else if (red<median) {
			red=red-detensifier;
		}
		if (green>median) {
			green=green+intensifier;
		} else if (red<median) {
			green=green-detensifier;
		}
		if (blue>median) {
			blue=blue+intensifier;
		} else if (red<median) {
			blue=blue-detensifier;
		}
		if (red>1.0f) red=1.0f;
		if (green>1.0f) green=1.0f;
		if (blue>1.0f) blue=1.0f;
		if (red<0.0f) red=0.0f;
		if (green<0.0f) green=0.0f;
		if (blue<0.0f) blue=0.0f;
		
		if (red>0.667f&&green>0.667f&&blue>0.667f) {
			float max=General.min(values);

			if (red<=max) red=red-0.5f;
			if (green<=max) green=green-0.5f;
			if (blue<=max) blue=blue-0.5f;
			
		}
		
		return new Color(red, green, blue);
	}
}
