package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TCharDoubleHashMap;

public class PeptideUtils {
	public static String getSmartDecoy(String peptide, byte charge, HashSet<String> backgroundProteome, SearchParameters parameters) {
		FragmentationModel model=new FragmentationModel(peptide, parameters.getAAConstants());
		double[] primaryIons=model.getPrimaryIons(parameters.getFragType(), charge);
		
		String decoy=reverse(peptide, parameters);
		int attempts=0;
		int maxTries=10;
		while (attempts<maxTries) {
			if (backgroundProteome.contains(decoy)) {
				decoy=shuffle(decoy, parameters);
			} else {
				model=new FragmentationModel(decoy, parameters.getAAConstants());
				double[] decoyIons=model.getPrimaryIons(parameters.getFragType(), charge);
				int matches=0;
				for (double decoyFragment : decoyIons) {
					Optional<Double> match=parameters.getFragmentTolerance().getMatch(primaryIons, decoyFragment);
					if (match.isPresent()) matches++;
				}
				
				float fraction=matches/(float)primaryIons.length;
				if (fraction<=0.4f) {
					break;
				} else {
					// otherwise too much overlap
					decoy=shuffle(decoy, parameters);
				}
			}
			attempts++;
		}
		return decoy;
	}
	
	public static String getDecoy(String peptide, HashSet<String> backgroundProteome, SearchParameters parameters) {
		String decoy=reverse(peptide, parameters);
		int attempts=0;
		int maxTries=3;
		while (attempts<maxTries) {
			if (backgroundProteome.contains(decoy)) {
				decoy=shuffle(decoy, parameters);
			} else {
				break;
			}
			attempts++;
		}
		return decoy;
	}
	
	public static String reverse(String peptide, SearchParameters parameters) {
		Triplet<double[], double[], String[]>triplet=getMasses(peptide, parameters.getAAConstants());
		String[] aas=triplet.z;
		reverse(aas, parameters.getEnzyme());
		return getSequence(aas);
	}
	
	public static void reverse(String[] aas, DigestionEnzyme enzyme) {
		int start=0;
		if (enzyme.isTargetPostSite(aas[start].charAt(0))) {
			start++;
		}
		int stop=aas.length-1;
		if (enzyme.isTargetPreSite(aas[stop].charAt(0))) {
			stop--;
		}
		
		while (start<=stop) {
			String c=aas[start];
			aas[start]=aas[stop];
			aas[stop]=c;
			start++;
			stop--;
		}
	}

	/**
	 * generates reliable random shuffle
	 * @param peptide
	 * @param enzyme
	 * @return
	 */
	public static String shuffle(String peptide, SearchParameters parameters) {
		Triplet<double[], double[], String[]>triplet=getMasses(peptide, parameters.getAAConstants());
		String[] aas=triplet.z;
		shuffle(aas, 0, parameters.getEnzyme());
		return getSequence(aas);
	}

	/**
	 * generates reliable random shuffle
	 * @param peptide
	 * @param enzyme
	 * @return
	 */
	public static String shuffle(String peptide, int shuffleSeed, SearchParameters parameters) {
		Triplet<double[], double[], String[]>triplet=getMasses(peptide, parameters.getAAConstants());
		String[] aas=triplet.z;
		shuffle(aas, shuffleSeed, parameters.getEnzyme());
		return getSequence(aas);
	}
	
	public static void shuffle(String[] aas, int shuffleSeed, DigestionEnzyme enzyme) {
		int start=0;
		if (enzyme.isTargetPostSite(aas[start].charAt(0))) {
			start++;
		}
		int stop=aas.length-1;
		if (enzyme.isTargetPreSite(aas[stop].charAt(0))) {
			stop--;
		}
		int diff=(stop)-start;
		
		// String.hashCode() is cross-platform consistent
		int seed=RandomGenerator.randomInt(shuffleSeed)+getSequence(aas).hashCode();
		for (int i=0; i<aas.length; i++) {
			seed=RandomGenerator.randomInt(seed);
			int index1=start+Math.abs(seed%diff+1);
			
			seed=RandomGenerator.randomInt(seed);
			int index2=start+Math.abs(seed%diff+1);
			if (index1!=index2) {
				String c=aas[index1];
				aas[index1]=aas[index2];
				aas[index2]=c;
			}
		}
	}
	
	public static String getSequence(String[] aas) {
		StringBuilder sb=new StringBuilder();
		for (String aa : aas) {
			sb.append(aa);
		}
		return sb.toString();
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
				String aaString=aas.get(masses.size()-1);
				neutralLosses.set(masses.size()-1, MassConstants.getNeutralLoss(aaString.charAt(0), modificationMass));
				aas.set(masses.size()-1, aaString+(modificationMass>=0?"[+":"[")+modificationMass+"]");
			} else {
				masses.add(aaConstants.getMass(ca[i]));
				neutralLosses.add(0.0);
				aas.add(Character.toString(ca[i]));
			}
		}
		return new Triplet<double[], double[], String[]>(masses.toArray(), neutralLosses.toArray(), aas.toArray(new String[aas.size()]));
	}
	
	public static int getNumberOfMods(String sequence, int nominalMass) {
		char[] ca=sequence.toCharArray();

		int total=0;
		for (int i = 0; i < ca.length; i++) {
			if (ca[i]=='[') {
				StringBuilder sb=new StringBuilder();
				i++;
				while (ca[i]!=']') {
					sb.append(ca[i]);
					i++;
				}
				String massText = sb.toString();
				double modificationMass = Double.valueOf(massText);
				
				if (Math.round(modificationMass)==nominalMass) {
					total++;
				}
			}
		}
		return total;
	}

	private static final DecimalFormat SKYLINE_DF = new DecimalFormat(".#");
	private static final DecimalFormat SKYLINE_PEAK_BOUNDARIES_DF = new DecimalFormat("#");

	public static String formatForSkyline(String sequence, AminoAcidConstants aaConstants) {
		return formatForSkyline(sequence, aaConstants, SKYLINE_DF);
	}
	
	public static String formatForSkylinePeakBoundaries(String sequence, AminoAcidConstants aaConstants) {
		return formatForSkyline(sequence, aaConstants, SKYLINE_PEAK_BOUNDARIES_DF);
	}
	
	public static String formatForEncyclopeDIA(String sequence, AminoAcidConstants aaConstants) {
		return formatForSkyline(sequence, aaConstants, null);
	}
	
	public static String formatForSkyline(String sequence, AminoAcidConstants aaConstants, DecimalFormat df) {
		char[] ca=sequence.toCharArray();
		TCharDoubleHashMap fixedMods=aaConstants.getFixedMods();
		
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
				String formattedMass=df==null?massText:df.format(modificationMass);
				if (formattedMass.charAt(0)!='+'&&formattedMass.charAt(0)!='-') {
					formattedMass="+"+formattedMass;
				}
				aas.set(aas.size()-1, aas.get(aas.size()-1)+"["+formattedMass+"]");
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
					double mass=fixedMods.get(aaChar);
					if (mass!=0.0f) {
						String formattedMass=df==null?Double.toString(mass):df.format(mass);
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
	
}
