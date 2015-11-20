package edu.washington.gs.maccoss.encyclopedia.algorithms;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;

public class IsotopicDistributionCalculator {
	private static final float[] hDist=new float[] {0.999855f, 0.000145f};
	private static final float[] cDist=new float[] {0.98916f, 0.01084f};
	private static final float[] oDist=new float[] {0.99757601f, 0.000378998479f, 0.002044992f};
	private static final float[] nDist=new float[] {0.99633f, 0.00366f};
	private static final float[] sDist=new float[] {0.95021f, 0.00745f, 0.04221f, 0.0f, 0.00013f};
	private static final float[][] abundance=new float[][] {hDist, cDist, oDist, nDist, sDist};

	
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
	public static float[] getIsotopeDistribution(String sequence, AminoAcidConstants constants) {
		int[] base=new int[] {2, 0, 1, 0, 0}; 
		for (char c : sequence.toCharArray()) {
			int[] aaProportion=constants.getAminoAcidProportions(c);
			if (aaProportion!=null) {
				for (int i=0; i<base.length; i++) {
					base[i]+=aaProportion[i];
				}
			}
		}
		return getIsotopeDistribution(base);
	}
	
	/**
	 * this algorithm follows Kubinyi, Analytica Chimica Acta, 247 (1991) 107-119
	 * @param counts
	 * @return
	 */
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
