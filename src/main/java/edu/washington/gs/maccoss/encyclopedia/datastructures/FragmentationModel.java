package edu.washington.gs.maccoss.encyclopedia.datastructures;

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
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;

//@Immutable
public class FragmentationModel {
	private final double[] masses;
	private final double[] neutralLosses;
	private final String[] aas;
	
	public FragmentationModel(String modifiedSequence, AminoAcidConstants aaConstants) {
		Triplet<double[], double[], String[]> tuple=PeptideUtils.getMasses(modifiedSequence, aaConstants);
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
	
	public AnnotatedLibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params) {
		return getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params, 0.0);
	}

	public AnnotatedLibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params, double minimumMass) {
		FragmentIon[] ions=getPrimaryIonObjects(params.getFragType(), precursorCharge);
		TDoubleArrayList ionsList=new TDoubleArrayList();
		ArrayList<FragmentIon> annotationList=new ArrayList<FragmentIon>();
		for (int i=0; i<ions.length; i++) {
			if (ions[i].mass>=minimumMass) {
				ionsList.add(ions[i].mass);
				annotationList.add(ions[i]);
			}
		}
		double[] masses=ionsList.toArray();
		
		float[] unitIntensities=new float[masses.length];
		Arrays.fill(unitIntensities, 1.0f);
		
		float[] unitCorrelation=new float[masses.length];
		Arrays.fill(unitCorrelation, 1.0f);

		String sequence=getModifiedSequence(params.getAAConstants());
		double precursorMZ=getChargedMass(precursorCharge);

		return new AnnotatedLibraryEntry(filename, accessions, 1, precursorMZ, precursorCharge, sequence, 1, retentionTime, 0.0f, masses, unitIntensities, unitCorrelation, annotationList.toArray(new FragmentIon[annotationList.size()]));
	}
	
	public String[] getAas() {
		return aas;
	}
	public String toString() {
		StringBuilder sb=new StringBuilder();
		for (String aa : aas) {
			sb.append(aa);
		}
		return sb.toString();
	}
	
	public static Pair<Character, Double> parseAA(String aa) {
		char c=aa.charAt(0);
		if (aa.length()>1) {
			double mod=Double.parseDouble(aa.substring(aa.indexOf('[')+1, aa.indexOf(']')));
			return new Pair<Character, Double>(c, mod);
		}
		return new Pair<Character, Double>(c, null);
	}
	
	public String getModifiedSequence(AminoAcidConstants aaConstants) {
		TCharDoubleHashMap fixedMods=aaConstants.getFixedMods();
		
		StringBuilder sb=new StringBuilder();
		for (String aa : aas) {
			sb.append(aa);

			if (aa.length()==1) {
				char aaChar=aa.charAt(0);
				if (fixedMods.contains(aaChar)) {
					double mass=fixedMods.get(aaChar);
					if (mass!=0.0f) {
						String formattedMass=Double.toString(mass);
						if (formattedMass.charAt(0)!='+'&&formattedMass.charAt(0)!='-') {
							formattedMass="+"+formattedMass;
						}
						sb.append("["+formattedMass+"]");
					}
				}
			}
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
			bs[i]=new FragmentIon(bs[i].mass+MassConstants.nh3, bs[i].index, IonType.c);
		}
		return bs;
	}
	
	public FragmentIon[] getZIons() {
		FragmentIon[] ys=getYIons();
		for (int i=0; i<ys.length; i++) {
			ys[i]=new FragmentIon(ys[i].mass-MassConstants.nh3, ys[i].index, IonType.z);
		}
		return ys;
	}
	
	public FragmentIon[] getZp1Ions() {
		FragmentIon[] ys=getYIons();
		for (int i=0; i<ys.length; i++) {
			ys[i]=new FragmentIon(ys[i].mass-MassConstants.nh3+MassConstants.hydrogenMass, ys[i].index, IonType.z1);
		}
		return ys;
	}

	public FragmentIon[] getBIons() {
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		
		ArrayList<FragmentIon> rolling=new ArrayList<FragmentIon>(); // seeds
		rolling.add(new FragmentIon(MassConstants.protonMass, (byte)0, IonType.b));
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
		rolling.add(new FragmentIon(MassConstants.oh2+MassConstants.protonMass, (byte)0, IonType.y));
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

}
