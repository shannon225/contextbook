package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.list.array.TDoubleArrayList;

//@Immutable
public class FragmentationModel {
	private final double[] masses;
	private final double[] neutralLosses;
	private final String[] aas;
	
	public FragmentationModel(String modifiedSequence, AminoAcidConstants aaConstants) {
		Triplet<double[], double[], String[]> tuple = getMasses(modifiedSequence, aaConstants);
		masses=tuple.x;
		neutralLosses=tuple.y;
		aas=tuple.z;
	}
	
	public LibraryEntry getUnitSpectrum(byte precursorCharge, SearchParameters params) {
		double[] ions=getPrimaryIons(params.getFragType(), precursorCharge);
		float[] unitIntensities=new float[ions.length];
		Arrays.fill(unitIntensities, 1.0f);

		String sequence=getModifiedSequence();
		double precursorMZ=params.getAAConstants().getChargedMass(sequence, precursorCharge);

		return new LibraryEntry(precursorMZ, precursorCharge, sequence, 1, 0.0f, 0.0f, ions, unitIntensities);
	}
	
	public String[] getAas() {
		return aas;
	}
	
	public static Pair<Character, Double> parseAA(String aa) {
		char c=aa.charAt(0);
		if (aa.length()>1) {
			double mod=Double.parseDouble(aa.substring(aa.indexOf('[')+1, aa.indexOf(']')));
			return new Pair<Character, Double>(c, mod);
		}
		return new Pair<Character, Double>(c, null);
	}
	
	public String getModifiedSequence() {
		StringBuilder sb=new StringBuilder();
		for (String string : aas) {
			sb.append(string);
		}
		return sb.toString();
	}
	
	/**
	 * returns sorted array of sprimary ions
	 * @param type
	 * @return
	 */
	public double[] getPrimaryIons(FragmentationType type, byte precursorCharge) {
		switch (type) {
			case YONLY:
				double[] yIons=getYIons();
				if (precursorCharge>2) {
					return concatAndSort(yIons, getPlus2s(yIons));
				} else {
					return yIons;
				}
			case CID:
				double[] yIonsCID=getYIons();
				double[] bIonsCID=getBIons();
				if (precursorCharge>2) {
					return concatAndSort(yIonsCID, getPlus2s(yIonsCID), bIonsCID, getPlus2s(bIonsCID));
				} else {
					return concatAndSort(bIonsCID, yIonsCID);
				}
			case ETD:
				double[] cIonsCID=getCIons();
				double[] zIonsCID=getZIons();
				double[] zp1IonsCID=getZp1Ions();
				if (precursorCharge>3) { // one charge gets quenched in fragmentation
					return concatAndSort(cIonsCID, getPlus2s(cIonsCID), zIonsCID, getPlus2s(zIonsCID), zp1IonsCID, getPlus2s(zp1IonsCID));
				} else {
					return concatAndSort(getCIons(), getZIons(), getZp1Ions());
				}
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+type+"]");
		}
	}
	
	public static double[] getPlus2s(double[] masses) {
		double[] p2=new double[masses.length];
		for (int i=0; i<p2.length; i++) {
			p2[i]=(masses[i]+MassConstants.protonMass)/2.0;
		}
		return p2;
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
	public static Triplet<double[], double[], String[]> getMasses(String sequence, AminoAcidConstants aaConstants) {
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
					masses.add(aaConstants.getMass(ca[i]));
					neutralLosses.add(0.0);
					aas.add(Character.toString(ca[i]));
				}
				String massText = sb.toString();
				double modificationMass = Double.valueOf(massText);
				masses.set(masses.size()-1, masses.get(masses.size()-1)+modificationMass);
				neutralLosses.set(masses.size()-1, MassConstants.getNeutralLoss(modificationMass));
				aas.set(masses.size()-1, aas.get(masses.size()-1)+(modificationMass>=0?"[+":"[")+modificationMass+"]");
			} else {
				masses.add(aaConstants.getMass(ca[i]));
				neutralLosses.add(0.0);
				aas.add(Character.toString(ca[i]));
			}
		}
		return new Triplet<double[], double[], String[]>(masses.toArray(), neutralLosses.toArray(), aas.toArray(new String[aas.size()]));
	}
	
}
