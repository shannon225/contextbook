package edu.washington.gs.maccoss.encyclopedia.algorithms;

import gnu.trove.map.hash.TCharObjectHashMap;

public class IsotopicDistributionCalculator {
	private static final TCharObjectHashMap<int[]> map=new TCharObjectHashMap<int[]>();
	private static final float[] hDist=new float[] {0.999855f, 0.000145f};
	private static final float[] cDist=new float[] {0.98916f, 0.01084f};
	private static final float[] oDist=new float[] {0.99757601f, 0.000378998479f, 0.002044992f};
	private static final float[] nDist=new float[] {0.99633f, 0.00366f};
	private static final float[] sDist=new float[] {0.95021f, 0.00745f, 0.04221f, 0.0f, 0.00013f};
	private static final float[][] abundance=new float[][] {hDist, cDist, oDist, nDist, sDist};

	static {
		map.put('A', new int[] {5, 3, 1, 1, 0});
		map.put('C', new int[] {8, 5, 2, 2, 1}); // assumes +57 alkylation
		map.put('D', new int[] {5, 4, 3, 1, 0});
		map.put('E', new int[] {7, 5, 3, 1, 0});
		map.put('F', new int[] {9, 9, 1, 1, 0});
		map.put('G', new int[] {3, 2, 1, 1, 0});
		map.put('H', new int[] {7, 6, 1, 3, 0});
		map.put('I', new int[] {11, 6, 1, 1, 0});
		map.put('K', new int[] {12, 6, 1, 2, 0});
		map.put('L', new int[] {11, 6, 1, 1, 0});
		map.put('M', new int[] {9, 5, 1, 1, 1});
		map.put('N', new int[] {6, 4, 2, 2, 0});
		map.put('P', new int[] {7, 5, 1, 1, 0});
		map.put('Q', new int[] {8, 5, 2, 2, 0});
		map.put('R', new int[] {12, 6, 1, 4, 0});
		map.put('S', new int[] {5, 3, 2, 1, 0});
		map.put('T', new int[] {7, 4, 2, 1, 0});
		map.put('V', new int[] {9, 5, 1, 1, 0});
		map.put('W', new int[] {10, 11, 1, 2, 0});
		map.put('Y', new int[] {9, 9, 2, 1, 0});
	}
	
	public static float[] normalizeToMax(float[] values) {
		float max=0.0f;
		for (int i = 0; i < values.length; i++) {
			if (values[i]>max) {
				max=values[i];
			}
		}
		if (max==0.0f) return values;
		float[] ret=new float[values.length];
		for (int i = 0; i < values.length; i++) {
			ret[i]=values[i]/max;
		}
		return ret;
	}

	/**
	 * I ignore modifications in isotope distribution calculation. This is an oversight.
	 * @param sequence
	 * @return
	 */
	public static float[] getIsotopeDistribution(String sequence) {
		int[] base=new int[] {2, 0, 1, 0, 0}; 
		for (char c : sequence.toCharArray()) {
			int[] aaProportion=map.get(c);
			if (aaProportion!=null) {
				for (int i=0; i<base.length; i++) {
					base[i]+=aaProportion[i];
				}
			}
		}
		return getIsotopeDistribution(base);
	}
	
	static float[] getIsotopeDistribution(int[] counts) {
		int p=0;
		int q=0;
		
		float[] cpattern=new float[5];
		cpattern[0]=1.0f;
		
		for (int j=0; j<abundance.length; j++) { // types of atoms (1...j...n)
			for (int i=0; i<counts[j]; i++) { // number of atom(j) in molecule
				float[] dist=new float[cpattern.length];
				
				for (int k=p; k<=q; k++) { // calculate the isotope distribution
					for (int l=0; l<abundance[j].length; l++) { // isotopes of atom(j)
						if (k+l>=cpattern.length) break;
						dist[k+l]+=cpattern[k]*abundance[j][l];
					}
				}
				
				q=q+abundance[j].length-1;
				if (q>=cpattern.length) {
					q=cpattern.length-1;
				}
				
				float max=0;
				for (int k=p; k<=q; k++) { // calculate the isotope distribution
					if (dist[k]>max) {
						max=dist[k];
					}
				}
				
				if (max>0) {
					for (int k=p; k<=q; k++) { // calculate the isotope distribution
						dist[k]=dist[k]/max;
					}
				}

				// create new isotope pattern
				for (int k=p; k<=q; k++) { // calculate the isotope distribution
					cpattern[k]=dist[k];
				}
			}
		}
		
		return cpattern;
	}
}
