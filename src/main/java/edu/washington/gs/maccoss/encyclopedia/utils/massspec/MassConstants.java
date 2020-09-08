package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public class MassConstants {
	public final static double neutronMass=1.0086649158849;
	public final static double protonMass=1.00727646681290;
	public final static double hydrogenMass=1.007825032071;
	public static final double carbonMass=12.0000000000000;
	public static final double oxygenMass=15.9949146195616;
	public static final double nitrogenMass=14.00307400486;
	public static final double siliconMass=27.97692653505;
	public static final double sodiumMass=22.989769282019;
	public final static double oh2=oxygenMass+2*hydrogenMass;
	public final static double nh3=nitrogenMass+3*hydrogenMass;
	public final static double co=carbonMass+oxygenMass;

	public static double getPeptideMass(double chargedMass, byte charge) {
		return chargedMass*charge-protonMass*charge;
	}

	public static double getChargedIsotopeMass(double precursorMz, byte charge, byte isotope) {
		return precursorMz + (isotope * neutronMass / charge);
	}
}
