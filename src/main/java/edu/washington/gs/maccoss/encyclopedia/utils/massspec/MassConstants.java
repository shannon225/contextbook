package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PeptideModification;
import edu.washington.gs.maccoss.encyclopedia.datastructures.AminoAcidConstants;

public class MassConstants {
	public final static double neutronMass=1.0086649158849;
	public final static double protonMass=1.00727646681290;
	public final static double hydrogenMass=1.007825032071;
	public static final double carbonMass=12.0000000000000;
	public static final double oxygenMass=15.9949146195616;
	public static final double nitrogenMass=14.00307400486;
	public final static double oh2=oxygenMass+2*hydrogenMass;
	public final static double nh3=nitrogenMass+3*hydrogenMass;
	public final static double co=carbonMass+oxygenMass;
	
	public static double getAccurateModificationMass(char aa, double modificationMass) {
		if (aa=='C') {
			if (tolerance.equals(57.0, modificationMass)) { // Carbamidomethyl
				return 57.0214635;
			} else if (tolerance.equals(58.0, modificationMass)) { // Carboxymethyl
					return 58.005479;
			} else if (tolerance.equals(46.0, modificationMass)) { // MMTS
				return 45.987721;
			} else if (tolerance.equals(99.0, modificationMass)) { // Carbamidomethyl + acetyl
				return 57.0214635+42.010565;
			} else if (tolerance.equals(40.0, modificationMass)) { // Carbamidomethyl - pyro-glu
				return 57.0214635-17.026549;
			}
		}
		
		if (aa=='M'||aa=='W') {
			if (tolerance.equals(16.0, modificationMass)) { // Oxidation
				return 15.994915;
			} else if (tolerance.equals(58.0, modificationMass)) { // Ox + acetyl
				return 42.010565+15.994915;
			}
		}
		
		if (aa=='Q') {
			if (tolerance.equals(-17.0, modificationMass)) { // pyro-glu
				return -17.026549;
			}
		}
		
		if (aa=='S'||aa=='T'||aa=='Y') {
			if (tolerance.equals(80.0, modificationMass)) { // Phospho
				return 79.966331;
			} else if (tolerance.equals(122.0, modificationMass)) { // Phospho + acetyl
				return 42.010565+79.966331;
			}
		}

		if (tolerance.equals(42.0, modificationMass)) { // acetyl
			return 42.010565;
		}
		
		for (PeptideModification mod : AminoAcidConstants.getDefaultLocalizationModifications()) {
			if (mod.isModificationMass(aa, modificationMass)) {
				return mod.getMass();
			}
		}
		
		return modificationMass;
	}
	
	private static final MassTolerance tolerance=new MassTolerance(1.0); // 1 ppm is about the accuracy of floats 
	
	public static double getPeptideMass(double chargedMass, byte charge) {
		return chargedMass*charge-protonMass*charge;
	}

	public static double getChargedIsotopeMass(double precursorMz, byte charge, byte isotope) {
		return precursorMz + (isotope * neutronMass / charge);
	}
}
