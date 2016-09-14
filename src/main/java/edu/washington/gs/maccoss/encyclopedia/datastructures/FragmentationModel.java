package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TCharFloatHashMap;

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

	public double getChargedMass(byte charge) {
		double mass=MassConstants.oh2;
		for (int i=0; i<masses.length; i++) {
			mass+=masses[i];
		}
		return (mass+MassConstants.protonMass*charge)/charge;
	}
	
	public LibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params) {
		return getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params, 0.0);
	}

	public AnnotatedLibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params, double minimumMass) {
		FragmentIon[] ions=getPrimaryIonObjects(params.getFragType(), precursorCharge);
		TDoubleArrayList ionsList=new TDoubleArrayList();
		ArrayList<String> annotationList=new ArrayList<String>();
		for (int i=0; i<ions.length; i++) {
			if (ions[i].mass>=minimumMass) {
				ionsList.add(ions[i].mass);
				annotationList.add(ions[i].toString());
			}
		}
		double[] masses=ionsList.toArray();
		
		float[] unitIntensities=new float[masses.length];
		Arrays.fill(unitIntensities, 1.0f);
		
		float[] unitCorrelation=new float[masses.length];
		Arrays.fill(unitCorrelation, 1.0f);

		String sequence=getModifiedSequence();
		double precursorMZ=getChargedMass(precursorCharge);

		return new AnnotatedLibraryEntry(filename, accessions, 1, precursorMZ, precursorCharge, sequence, 1, retentionTime, 0.0f, masses, unitIntensities, unitCorrelation, annotationList.toArray(new String[annotationList.size()]));
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
		FragmentIon[] ions=getPrimaryIonObjects(type, precursorCharge);
		double[] masses=new double[ions.length];
		for (int i=0; i<ions.length; i++) {
			masses[i]=ions[i].mass;
		}
		return masses;
	}

	public FragmentIon[] getPrimaryIonObjects(FragmentationType type, byte precursorCharge) {
		switch (type) {
			case YONLY:
				FragmentIon[] yIons=getYIons();
				if (precursorCharge>2) {
					return concatAndSort(yIons, getPlus2s(yIons));
				} else {
					return yIons;
				}
			case CID:
				FragmentIon[] yIonsCID=getYIons();
				FragmentIon[] bIonsCID=getBIons();
				if (precursorCharge>2) {
					return concatAndSort(yIonsCID, getPlus2s(yIonsCID), bIonsCID, getPlus2s(bIonsCID));
				} else {
					return concatAndSort(bIonsCID, yIonsCID);
				}
			case ETD:
				FragmentIon[] cIonsCID=getCIons();
				FragmentIon[] zIonsCID=getZIons();
				FragmentIon[] zp1IonsCID=getZp1Ions();
				if (precursorCharge>3) { // one charge gets quenched in fragmentation
					return concatAndSort(cIonsCID, getPlus2s(cIonsCID), zIonsCID, getPlus2s(zIonsCID), zp1IonsCID, getPlus2s(zp1IonsCID));
				} else {
					return concatAndSort(cIonsCID, zIonsCID, zp1IonsCID);
				}
			default:
				throw new EncyclopediaException("Unknown fragmentation type ["+type+"]");
		}
	}
	
	public static FragmentIon[] getPlus2s(FragmentIon[] masses) {
		FragmentIon[] p2=new FragmentIon[masses.length];
		for (int i=0; i<p2.length; i++) {
			p2[i]=new FragmentIon((masses[i].mass+MassConstants.protonMass)/2.0, masses[i].index, IonType.getPlus2(masses[i].type));
		}
		return p2;
	}

	private static FragmentIon[] concatAndSort(FragmentIon[]... a) {
		int length=0;
		for (FragmentIon[] ds : a) {
			length+=ds.length;
		}
		FragmentIon[] c=new FragmentIon[length];
		int current=0;
		for (FragmentIon[] ds : a) {
			System.arraycopy(ds, 0, c, current, ds.length);
			current+=ds.length;
		}
		Arrays.sort(c);
		return c;
	}

	public FragmentIon[] getCIons() {
		FragmentIon[] bs=getBIons();
		for (int i=0; i<bs.length; i++) {
			bs[i]=new FragmentIon(bs[i].mass+17.02654911, bs[i].index, IonType.c);
		}
		return bs;
	}
	
	public FragmentIon[] getZIons() {
		FragmentIon[] ys=getYIons();
		for (int i=0; i<ys.length; i++) {
			ys[i]=new FragmentIon(ys[i].mass-17.02654911, ys[i].index, IonType.z);
		}
		return ys;
	}
	
	public FragmentIon[] getZp1Ions() {
		FragmentIon[] ys=getYIons();
		for (int i=0; i<ys.length; i++) {
			ys[i]=new FragmentIon(ys[i].mass-16.01872407, ys[i].index, IonType.z1);
		}
		return ys;
	}

	public FragmentIon[] getBIons() {
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		
		ArrayList<FragmentIon> rolling=new ArrayList<FragmentIon>(); // seeds
		rolling.add(new FragmentIon(1.007276467, (byte)0, IonType.b));
		for (byte i = 0; i < masses.length; i++) {
			int index=i;
			ArrayList<FragmentIon> neutrals=new ArrayList<FragmentIon>();
			for (int j = 0; j < rolling.size(); j++) {
				rolling.set(j, rolling.get(j).increment(masses[index]));
				ions.add(rolling.get(j));
				if (neutralLosses[index]>0.0) {
					FragmentIon nl=rolling.get(j).neutralLoss(neutralLosses[index]);
					neutrals.add(nl);
					ions.add(nl);
				}
			}
			if (neutrals.size()>0) {
				rolling.addAll(neutrals);
			}
		}
		FragmentIon[] ionArray = ions.toArray(new FragmentIon[ions.size()]);
		Arrays.sort(ionArray);
		return ionArray;
	}
	
	public FragmentIon[] getYIons() {
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		
		ArrayList<FragmentIon> rolling=new ArrayList<FragmentIon>();
		rolling.add(new FragmentIon(19.01784117, (byte)0, IonType.y));
		for (int i = 0; i < masses.length; i++) {
			int index=masses.length-1-i;
			ArrayList<FragmentIon> neutrals=new ArrayList<FragmentIon>();
			for (int j = 0; j < rolling.size(); j++) {
				rolling.set(j, rolling.get(j).increment(masses[index]));
				ions.add(rolling.get(j));
				if (neutralLosses[index]>0.0) {
					FragmentIon nl=rolling.get(j).neutralLoss(neutralLosses[index]);
					neutrals.add(nl);
					ions.add(nl);
				}
			}
			if (neutrals.size()>0) {
				rolling.addAll(neutrals);
			}
		}
		FragmentIon[] ionArray = ions.toArray(new FragmentIon[ions.size()]);
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

	private static final DecimalFormat SKYLINE_DF = new DecimalFormat(".#");
	private static final DecimalFormat SKYLINE_PEAK_BOUNDARIES_DF = new DecimalFormat("#");

	public static String formatForSkyline(String sequence, AminoAcidConstants aaConstants) {
		return formatForSkyline(sequence, aaConstants, SKYLINE_DF);
	}
	
	public static String formatForSkylinePeakBoundaries(String sequence, AminoAcidConstants aaConstants) {
		return formatForSkyline(sequence, aaConstants, SKYLINE_PEAK_BOUNDARIES_DF);
	}
	
	public static String formatForSkyline(String sequence, AminoAcidConstants aaConstants, DecimalFormat df) {
		char[] ca=sequence.toCharArray();
		TCharFloatHashMap fixedMods=aaConstants.getFixedMods();
		
		ArrayList<String> aas=new ArrayList<String>();
		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='[') {
				StringBuilder sb=new StringBuilder();
				i++;
				while (ca[i]!=']') {
					sb.append(ca[i]);
					i++;
				}
				if (aas.size()==0) {
					// handling of n-termini mods assumes you can't have multiple []s in a row
					i++;
					aas.add(Character.toString(ca[i]));
				}
				String massText = sb.toString();
				double modificationMass = Double.valueOf(massText);
				aas.set(aas.size()-1, aas.get(aas.size()-1)+(modificationMass>=0?"[+":"[")+df.format(modificationMass)+"]");
			} else {
				aas.add(Character.toString(ca[i]));
			}
		}
		
		StringBuilder sb=new StringBuilder();
		for (String aa : aas) {
			sb.append(aa);
			if (aa.length()==1) {
				char aaChar=aa.charAt(0);
				if (fixedMods.contains(aaChar)) {
					float mass=fixedMods.get(aaChar);
					if (mass!=0.0f) {
						sb.append((mass>=0?"[+":"[")+df.format(mass)+"]");
					}
				}
			}
		}
		return sb.toString();
	}
	
}
