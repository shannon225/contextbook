package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizer;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentationType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.IonType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassConstants;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import gnu.trove.list.array.TDoubleArrayList;

//@Immutable
public class FragmentationModel {
	private static final int NUMBER_OF_NEUTRONS_TO_CONSIDER_STILL_IN_RANGE=4;
	private final double[] masses;
	private final double[] modificationMasses;
	private final double[] neutralLosses;
	private final String[] aas;
	
	public FragmentationModel(double[] masses, double[] modificationMasses, double[] neutralLosses, String[] aas) {
		this.masses=masses;
		this.modificationMasses=modificationMasses;
		this.neutralLosses=neutralLosses;
		this.aas=aas;
	}

	public static AnnotatedLibraryEntry generateEntry(String peptideModSeq, String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, boolean isDecoy, SearchParameters params) {
		FragmentationModel model=PeptideUtils.getPeptideModel(peptideModSeq, params.getAAConstants());
		return model.getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params, isDecoy);
	}

	public double getChargedMass(byte charge) {
		double mass=MassConstants.oh2;
		for (int i=0; i<masses.length; i++) {
			mass+=masses[i];
		}
		return (mass+MassConstants.protonMass*charge)/charge;
	}
	
	public AnnotatedLibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params) {
		return getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params, 0.0, false);
	} 
	public AnnotatedLibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params, boolean isDecoy) {
		return getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params, 0.0, isDecoy);
	}

	public AnnotatedLibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params, double minimumMass, boolean forQuant) {
		return getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params, minimumMass, false, forQuant);
	} 
	public AnnotatedLibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params, double minimumMass, boolean isDecoy, boolean forQuant) {
		return getUnitSpectrum(filename, accessions, precursorCharge, retentionTime, params, null, minimumMass, false, forQuant);
	}
	public AnnotatedLibraryEntry getUnitSpectrum(String filename, HashSet<String> accessions, byte precursorCharge, float retentionTime, SearchParameters params, double[] targetMasses, double minimumMass, boolean isDecoy, boolean forQuant) {
		String sequence=getModifiedSequence();
		double precursorMZ=getChargedMass(precursorCharge);
		FragmentIon[] ions=getPrimaryIonObjects(params.getFragType(), precursorCharge, forQuant);
		MassTolerance fragmentTolerance=params.getFragmentTolerance();
		ions = FragmentIon.getUniqueFragments(ions, fragmentTolerance);
		
		if (targetMasses!=null) {
			targetMasses=targetMasses.clone();
			Arrays.sort(targetMasses);
		}

		TDoubleArrayList ionsList=new TDoubleArrayList();
		ArrayList<FragmentIon> annotationList=new ArrayList<FragmentIon>();
		for (int i=0; i<ions.length; i++) {
			if (ions[i].mass>=minimumMass) {
				if (targetMasses==null||fragmentTolerance.getIndex(targetMasses, ions[i].mass).isPresent()) {
					ionsList.add(ions[i].mass);
					annotationList.add(ions[i]);
				}
			}
		}
		double[] masses=ionsList.toArray();
		
		float[] unitIntensities=new float[masses.length];
		Arrays.fill(unitIntensities, 1.0f);
		
		float[] unitCorrelation=new float[masses.length];
		Arrays.fill(unitCorrelation, 1.0f);

		return new AnnotatedLibraryEntry(filename, accessions, 1, precursorMZ, precursorCharge, sequence, 1, retentionTime, 0.0f, masses, unitIntensities, unitCorrelation, annotationList.toArray(new FragmentIon[annotationList.size()]), isDecoy, params.getAAConstants());
	}
	public double[] getMasses() {
		return masses;
	}
	
	public double[] getModificationMasses() {
		return modificationMasses;
	}
	
	public double[] getNeutralLosses() {
		return neutralLosses;
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
	
	public String getModifiedSequence() {
		StringBuilder sb=new StringBuilder();
		for (String aa : aas) {
			sb.append(aa);
		}
		return sb.toString();
	}
	
	/**
	 * returns sorted array of sprimary ions
	 * @param type
	 * @return
	 */
	public double[] getPrimaryIons(FragmentationType type, byte precursorCharge, boolean forQuant) {
		FragmentIon[] ions=getPrimaryIonObjects(type, precursorCharge, forQuant);
		double[] masses=new double[ions.length];
		for (int i=0; i<ions.length; i++) {
			masses[i]=ions[i].mass;
		}
		return masses;
	}
	
	/**
	 * finds the ions that uniquely describe this model if it is modified in a
	 * way where ions might appear in the same precursor isolation window. If it
	 * is not, then returns Optional.empty() and you don't need to worry about
	 * modifications.
	 * 
	 * @param precursorRange
	 * @param type
	 * @param precursorCharge
	 * @param forQuant
	 * @return
	 */
	public Optional<FragmentIon[]> getModificationSpecificIonObjects(Range precursorRange, FragmentationType type, byte precursorCharge, boolean forQuant) {
		HashMap<String, FragmentationModel> availableModels=new HashMap<>();
		
		double precursorMz=getChargedMass(precursorCharge);
		for (int i=0; i<modificationMasses.length; i++) {
			double unmodifiedMass=precursorMz-modificationMasses[i]/precursorCharge;
			Range unmodifiedPrecursorRange=new Range((float)(unmodifiedMass-(NUMBER_OF_NEUTRONS_TO_CONSIDER_STILL_IN_RANGE*MassConstants.neutronMass/precursorCharge)), (float)unmodifiedMass);
			
			if (modificationMasses[i]!=0.0&&precursorRange.contains(unmodifiedPrecursorRange)) {
				String[] altAAs=aas.clone();
				double[] altMasses=masses.clone();
				double[] altModMasses=modificationMasses.clone();
				double[] altNLs=neutralLosses.clone();
				
				altAAs[i]=aas[i].substring(0, 1);
				altMasses[i]=masses[i]-modificationMasses[i];
				altModMasses[i]=0.0;
				altNLs[i]=0.0;
				
				FragmentationModel altModel=new FragmentationModel(altMasses, altModMasses, altNLs, altAAs);
				availableModels.put(altModel.getModifiedSequence(), altModel);
			}
		}
		
		if (availableModels.size()==0) return Optional.empty();
		
		final String peptideModSeq=getModifiedSequence();
		availableModels.put(peptideModSeq, this);
		return Optional.of(PhosphoLocalizer.getUniqueFragmentIons(peptideModSeq, precursorCharge, availableModels, type));
	}

	public FragmentIon[] getPrimaryIonObjects(FragmentationType type, byte precursorCharge, boolean forQuant) {
		return getPrimaryIonObjects(type, precursorCharge, true, forQuant);
	}
	public FragmentIon[] getPrimaryIonObjects(FragmentationType type, byte precursorCharge, boolean useNeutralLosses, boolean forQuant) {
		switch (type) {
			case HCD:
				if (forQuant) {
					// include B ions too
					FragmentIon[] yIonsCID=getYIons(useNeutralLosses);
					FragmentIon[] bIonsCID=getBIons(useNeutralLosses);
					if (precursorCharge>2) {
						return concatAndSort(yIonsCID, getPlus2s(yIonsCID), bIonsCID, getPlus2s(bIonsCID));
					} else {
						return concatAndSort(bIonsCID, yIonsCID);
					}
				} else {
					FragmentIon[] yIons=getYIons(useNeutralLosses);
					if (precursorCharge>2) {
						return concatAndSort(yIons, getPlus2s(yIons));
					} else {
						return yIons;
					}
				}
			case CID:
				FragmentIon[] yIonsCID=getYIons(useNeutralLosses);
				FragmentIon[] bIonsCID=getBIons(useNeutralLosses);
				if (precursorCharge>2) {
					return concatAndSort(yIonsCID, getPlus2s(yIonsCID), bIonsCID, getPlus2s(bIonsCID));
				} else {
					return concatAndSort(bIonsCID, yIonsCID);
				}
			case ETD:
				FragmentIon[] cIonsCID=getCIons(useNeutralLosses);
				FragmentIon[] zIonsCID=getZIons(useNeutralLosses);
				FragmentIon[] zp1IonsCID=getZp1Ions(useNeutralLosses);
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
		return getCIons(true);
	}

	public FragmentIon[] getCIons(boolean useNeutralLosses) {
		FragmentIon[] bs=getBIons(useNeutralLosses);
		for (int i=0; i<bs.length; i++) {
			bs[i]=new FragmentIon(bs[i].mass+MassConstants.nh3, bs[i].index, IonType.c);
		}
		return bs;
	}
	
	public FragmentIon[] getZIons() {
		return getZIons(true);
	}

	public FragmentIon[] getZIons(boolean useNeutralLosses) {
		FragmentIon[] ys=getYIons(useNeutralLosses);
		for (int i=0; i<ys.length; i++) {
			ys[i]=new FragmentIon(ys[i].mass-MassConstants.nh3, ys[i].index, IonType.z);
		}
		return ys;
	}
	
	public FragmentIon[] getZp1Ions() {
		return getZp1Ions(true);
	}

	public FragmentIon[] getZp1Ions(boolean useNeutralLosses) {
		FragmentIon[] ys=getYIons(useNeutralLosses);
		for (int i=0; i<ys.length; i++) {
			ys[i]=new FragmentIon(ys[i].mass-MassConstants.nh3+MassConstants.hydrogenMass, ys[i].index, IonType.z1);
		}
		return ys;
	}

	public FragmentIon[] getBIons() {
		return getBIons(true);
	}

	public FragmentIon[] getBIons(boolean useNeutralLosses) {
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		
		ArrayList<FragmentIon> rolling=new ArrayList<FragmentIon>(); // seeds
		rolling.add(new FragmentIon(MassConstants.protonMass, (byte)0, IonType.b));
		for (byte i = 0; i < masses.length; i++) {
			int index=i;
			ArrayList<FragmentIon> neutrals=new ArrayList<FragmentIon>();
			for (int j = 0; j < rolling.size(); j++) {
				rolling.set(j, rolling.get(j).increment(masses[index]));
				ions.add(rolling.get(j));
				if (useNeutralLosses&&neutralLosses[index]>0.0) {
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
		return getYIons(true);
	}

	public FragmentIon[] getYIons(boolean useNeutralLosses) {
		ArrayList<FragmentIon> ions=new ArrayList<FragmentIon>();
		
		ArrayList<FragmentIon> rolling=new ArrayList<FragmentIon>();
		rolling.add(new FragmentIon(MassConstants.oh2+MassConstants.protonMass, (byte)0, IonType.y));
		for (int i = 0; i < masses.length; i++) {
			int index=masses.length-1-i;
			ArrayList<FragmentIon> neutrals=new ArrayList<FragmentIon>();
			for (int j = 0; j < rolling.size(); j++) {
				rolling.set(j, rolling.get(j).increment(masses[index]));
				ions.add(rolling.get(j));
				if (useNeutralLosses&&neutralLosses[index]>0.0) {
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
