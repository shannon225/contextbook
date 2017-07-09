package edu.washington.gs.maccoss.encyclopedia.utils.math;

import java.util.Arrays;

import gnu.trove.map.hash.TDoubleDoubleHashMap;

public class BenjaminiHochberg {
	public static double[] calculateAdjustedPValues(double[] pValues) {
		double[] copy=pValues.clone();
		double[] adjusted=new double[copy.length];
		Arrays.sort(copy);
		
		TDoubleDoubleHashMap adjustedByP=new TDoubleDoubleHashMap();
		for (int i=copy.length-1; i>=0; i--) {
			if (i==copy.length-1) {
				adjusted[i]=copy[i];
			} else {
				double originalP=copy[i];
				double left=adjusted[i+1];
				double right=(copy.length/(double)(i+1))*originalP;
				adjusted[i]=Math.min(left, right);
				adjustedByP.put(originalP, Math.min(left, right));
			}
		}
		
		double[] ret=new double[pValues.length];
		for (int i=0; i<pValues.length; i++) {
			ret[i]=adjustedByP.get(pValues[i]);
		}
		return ret;
	}
}
