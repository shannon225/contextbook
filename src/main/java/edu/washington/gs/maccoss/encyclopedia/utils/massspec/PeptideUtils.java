package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

public class PeptideUtils {
	public static byte getExpectedChargeState(String peptide) {
		byte charge=1;
		for (int i=0; i<peptide.length(); i++) {
			switch (peptide.charAt(i)) {
				case 'K':
				case 'R':
				case 'H':
					charge++;
					break;
			}
		}
		return charge;
	}
	
	public static String getSmartDecoy(String peptide, byte charge, HashSet<String> backgroundProteome, SearchParameters parameters) {
		FragmentationModel model=PeptideUtils.getPeptideModel(peptide, parameters.getAAConstants());
		double[] primaryIons=model.getPrimaryIons(parameters.getFragType(), charge, false);
		
		String decoy=reverse(peptide, parameters.getAAConstants());
		int attempts=0;
		int maxTries=10;
		while (attempts<maxTries) {
			if (backgroundProteome.contains(decoy)) {
				decoy=shuffle(decoy, parameters);
			} else {
				model=PeptideUtils.getPeptideModel(decoy, parameters.getAAConstants());
				double[] decoyIons=model.getPrimaryIons(parameters.getFragType(), charge, false);
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
		String decoy=reverse(peptide, parameters.getAAConstants());
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

	public static String reverse(String peptide, AminoAcidConstants aminoAcidConstants) {
		String[] aas=getAAs(peptide, aminoAcidConstants);
		reverse(aas);
		return getSequence(aas);
	}
	
	/**
	 * keeps termini in place
	 * @param aas
	 */
	public static void reverse(String[] aas) {
		int start=1;
		int stop=aas.length-2;
		
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
		String[] aas=getPeptideModel(peptide, parameters.getAAConstants()).getAas();
		shuffle(aas, 0, parameters.getEnzyme());
		return getSequence(aas);
	}

	/**
	 * keeps termini in place, generates reliable random shuffle
	 * @param peptide
	 * @param enzyme
	 * @return
	 */
	public static String shuffle(String peptide, int shuffleSeed, SearchParameters parameters) {
		String[] aas=getPeptideModel(peptide, parameters.getAAConstants()).getAas();
		shuffle(aas, shuffleSeed, parameters.getEnzyme());
		return getSequence(aas);
	}
	
	public static void shuffle(String[] aas, int shuffleSeed, DigestionEnzyme enzyme) {
		int start=0;
		//if (enzyme.isTargetPostSite(aas[start].charAt(0))) {
			start++;
		//}
		int stop=aas.length-1;
		//if (enzyme.isTargetPreSite(aas[stop].charAt(0))) {
			stop--;
		//}
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

	public static String getCorrectedMasses(String sequence, AminoAcidConstants aminoAcidConstants) {
		char[] ca=sequence.toCharArray();
		
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
				String aaString=aas.get(aas.size()-1);
				char aaChar=aaString.charAt(0);
				modificationMass=aminoAcidConstants.getAccurateModificationMass(aaChar, modificationMass);

				aas.set(aas.size()-1, aaString+(modificationMass>=0?"[+":"[")+modificationMass+"]");
			} else {
				aas.add(Character.toString(ca[i]));
			}
		}
		StringBuilder sb=new StringBuilder();
		for (String aa : aas) {
			sb.append(aa);
		}
		return sb.toString();
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
	public static FragmentationModel getPeptideModel(String sequence, AminoAcidConstants aaConstants) {
		char[] ca=sequence.toCharArray();
		
		TDoubleArrayList masses=new TDoubleArrayList();
		TDoubleArrayList modifications=new TDoubleArrayList();
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
					modifications.add(0.0);
					aas.add(Character.toString(ca[i]));
				}
				String massText = sb.toString();
				double modificationMass = Double.valueOf(massText);
				String aaString=aas.get(masses.size()-1);
				char aaChar=aaString.charAt(0);
				modificationMass=aaConstants.getAccurateModificationMass(aaChar, modificationMass);

				masses.set(masses.size()-1, masses.get(masses.size()-1)+modificationMass);
				double neutralLoss = 0.0f;//aaConstants.getNeutralLoss(aaChar, modificationMass);
				neutralLosses.set(masses.size()-1, neutralLoss);
				modifications.set(masses.size()-1, modificationMass);
				aas.set(masses.size()-1, aaString+(modificationMass>=0?"[+":"[")+modificationMass+"]");
			} else {
				masses.add(aaConstants.getMass(ca[i]));
//				if (ca[i]=='D') {
//					neutralLosses.add(MassConstants.oh2);
//				} else if (ca[i]=='E') {
//					neutralLosses.add(MassConstants.oh2);
//				} else if (ca[i]=='N') {
//					neutralLosses.add(MassConstants.nh3);
//				} else if (ca[i]=='Q') {
//					neutralLosses.add(MassConstants.nh3);
//				} else {
//					neutralLosses.add(0.0);
//				}

				neutralLosses.add(0.0);
				modifications.add(0.0);
				aas.add(Character.toString(ca[i]));
			}
		}
		return new FragmentationModel(masses.toArray(), modifications.toArray(), neutralLosses.toArray(), aas.toArray(new String[aas.size()]));
	}

	public static String[] getAAs(String sequence, AminoAcidConstants aminoAcidConstants) {
		char[] ca=sequence.toCharArray();
		
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
				String aaString=aas.get(aas.size()-1);
				char aaChar=aaString.charAt(0);
				modificationMass=aminoAcidConstants.getAccurateModificationMass(aaChar, modificationMass);

				aas.set(aas.size()-1, aaString+(modificationMass>=0?"[+":"[")+modificationMass+"]");
			} else {
				aas.add(Character.toString(ca[i]));
			}
		}
		return aas.toArray(new String[aas.size()]);
	}

	
	public static String getPeptideSeq(String peptideModSeq) {
		StringBuilder sb=new StringBuilder();
		for (char c : peptideModSeq.toCharArray()) {
			if (Character.isLetter(c)) {
				sb.append(c);
			}
		}
		return sb.toString();
	}
	
	public static int[] getModIndicies(String sequence, int nominalMass) {
		char[] ca=sequence.toCharArray();

		TIntArrayList total=new TIntArrayList();
		int index=0;
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
					total.add(index);
				}
			} else {
				index++;
			}
		}
		return total.toArray();
	}
	
	public static int getNumberOfMods(String sequence, int nominalMass) {
		return getModIndicies(sequence, nominalMass).length;
	}

	private static final DecimalFormat SKYLINE_DF = new DecimalFormat(".#");
	private static final DecimalFormat SKYLINE_PEAK_BOUNDARIES_DF = new DecimalFormat("#");

	public static String formatForSkyline(String sequence) {
		return formatForSkyline(sequence, SKYLINE_DF, true);
	}
	
	public static String formatForSkylinePeakBoundaries(String sequence) {
		return formatForSkyline(sequence, SKYLINE_PEAK_BOUNDARIES_DF, true);
	}
	
	public static String formatForEncyclopeDIA(String sequence) {
		return formatForSkyline(sequence, null, true);
	}
	
	public static String formatForTPP(String sequence) {
		// doesn't support n&c termini
		return formatForSkyline(sequence, SKYLINE_PEAK_BOUNDARIES_DF, false);
	}
	
	private static String formatForSkyline(String sequence, DecimalFormat df, boolean addPlus) {
		char[] ca=sequence.toCharArray();
		
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
				if (addPlus&&formattedMass.charAt(0)!='+'&&formattedMass.charAt(0)!='-') {
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
		}
		return sb.toString();
	}
	
}
