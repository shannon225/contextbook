package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.list.array.TDoubleArrayList;

public class FragmentationModel {
	private final double[] masses;
	private final double[] neutralLosses;
	private final String[] aas;
	
	public FragmentationModel(String modifiedSequence) {
		Triplet<double[], double[], String[]> tuple = getMasses(modifiedSequence);
		masses=tuple.x;
		neutralLosses=tuple.y;
		aas=tuple.z;
	}
	
	public String[] getAas() {
		return aas;
	}
	
	/**
	 * returns sorted array of sprimary ions
	 * @param type
	 * @return
	 */
	public double[] getPrimaryIons(FragmentationType type) {
		switch (type) {
			case CID:
				return concatAndSort(getBIons(), getYIons());
			case ETD:
				return concatAndSort(getCIons(), getZIons(), getZp1Ions());
			default:
				throw new IllegalArgumentException("Unknown fragmentation type ["+type+"]");
		}
	}

	private static double[] concatAndSort(double[]... a) {
		int length=0;
		for (double[] ds : a) {
			length+=ds.length;
		}
		double[] c=new double[length];
		int current=0;
		for (double[] ds : a) {
			System.arraycopy(ds, 0, c, current, ds.length);
			current+=ds.length;
		}
		Arrays.sort(c);
		return c;
	}

	public double[] getCIons() {
		double[] bs=getBIons();
		for (int i=0; i<bs.length; i++) {
			bs[i]+=17.02654911;
		}
		return bs;
	}
	
	public double[] getZIons() {
		double[] ys=getYIons();
		for (int i=0; i<ys.length; i++) {
			ys[i]-=17.02654911;
		}
		return ys;
	}
	
	public double[] getZp1Ions() {
		double[] ys=getYIons();
		for (int i=0; i<ys.length; i++) {
			ys[i]-=16.01872407;
		}
		return ys;
	}

	public double[] getBIons() {
		TDoubleArrayList ions=new TDoubleArrayList();
		
		TDoubleArrayList rolling=new TDoubleArrayList();
		rolling.add(1.007276467);
		for (int i = 0; i < masses.length; i++) {
			int index=i;
			TDoubleArrayList neutrals=new TDoubleArrayList();
			for (int j = 0; j < rolling.size(); j++) {
				rolling.set(j, rolling.get(j)+masses[index]);
				ions.add(rolling.get(j));
				if (neutralLosses[index]>0.0) {
					neutrals.add(rolling.get(j)-neutralLosses[index]);
					ions.add(rolling.get(j)-neutralLosses[index]);
				}
			}
			if (neutrals.size()>0) {
				rolling.addAll(neutrals);
			}
		}
		double[] ionArray = ions.toArray();
		Arrays.sort(ionArray);
		return ionArray;
	}
	
	public double[] getYIons() {
		TDoubleArrayList ions=new TDoubleArrayList();
		
		TDoubleArrayList rolling=new TDoubleArrayList();
		rolling.add(19.01784117);
		for (int i = 0; i < masses.length; i++) {
			int index=masses.length-1-i;
			TDoubleArrayList neutrals=new TDoubleArrayList();
			for (int j = 0; j < rolling.size(); j++) {
				rolling.set(j, rolling.get(j)+masses[index]);
				ions.add(rolling.get(j));
				if (neutralLosses[index]>0.0) {
					neutrals.add(rolling.get(j)-neutralLosses[index]);
					ions.add(rolling.get(j)-neutralLosses[index]);
				}
			}
			if (neutrals.size()>0) {
				rolling.addAll(neutrals);
			}
		}
		double[] ionArray = ions.toArray();
		Arrays.sort(ionArray);
		return ionArray;
	}

	/**
	 * triplet is (masses, neutral losses, AAs)
	 * @param sequence
	 * @return
	 */
	public static Triplet<double[], double[], String[]> getMasses(String sequence) {
		char[] ca=sequence.toCharArray();
		
		TDoubleArrayList masses=new TDoubleArrayList();
		TDoubleArrayList neutralLosses=new TDoubleArrayList();
		ArrayList<String> aas=new ArrayList<String>();
		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='[') {
				StringBuilder sb=new StringBuilder();
				i++;
				while (ca[i]!=']') {
					sb.append(ca[i]);
					i++;
				}
				if (masses.size()==0) {
					// handling of n-termini mods assumes you can't have multiple []s in a row
					i++;
					masses.add(MassConstants.getMass(ca[i]));
					neutralLosses.add(0.0);
					aas.add(Character.toString(ca[i]));
				}
				String massText = sb.toString();
				double modificationMass = Double.valueOf(massText);
				masses.set(masses.size()-1, masses.get(masses.size()-1)+modificationMass);
				neutralLosses.set(masses.size()-1, MassConstants.getNeutralLoss(modificationMass));
				aas.set(masses.size()-1, aas.get(masses.size()-1)+"["+massText+"]");
			} else {
				masses.add(MassConstants.getMass(ca[i]));
				neutralLosses.add(0.0);
				aas.add(Character.toString(ca[i]));
			}
		}
		return new Triplet<double[], double[], String[]>(masses.toArray(), neutralLosses.toArray(), aas.toArray(new String[aas.size()]));
	}
}
