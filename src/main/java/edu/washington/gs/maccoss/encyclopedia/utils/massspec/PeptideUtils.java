package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentationModel;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import edu.washington.gs.maccoss.encyclopedia.utils.math.RandomGenerator;

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
		Triplet<double[], double[], String[]>triplet=FragmentationModel.getMasses(peptide, parameters.getAAConstants());
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
		Triplet<double[], double[], String[]>triplet=FragmentationModel.getMasses(peptide, parameters.getAAConstants());
		String[] aas=triplet.z;
		shuffle(aas, parameters.getEnzyme());
		return getSequence(aas);
	}
	
	public static void shuffle(String[] aas, DigestionEnzyme enzyme) {
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
		int seed=getSequence(aas).hashCode();
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
	
	private static String getSequence(String[] aas) {
		StringBuilder sb=new StringBuilder();
		for (String aa : aas) {
			sb.append(aa);
		}
		return sb.toString();
	}
}
