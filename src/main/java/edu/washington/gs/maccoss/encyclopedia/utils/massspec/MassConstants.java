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
	
	public static float getEstimateIMS(double chargedMass, byte charge) {
		// predicts 1/k0
		// discussed here: https://github.com/lazear/sage/pull/98
		double sqMzOverCharge=(chargedMass*chargedMass) / charge;
		
		return (float)(-1.660e+00+ 
		        (-3.798e-01 * Math.log1p(chargedMass)) + 
		        (-2.389e-04 * chargedMass) + 
		        (3.957e-01 * Math.log1p(sqMzOverCharge)) + 
		        (4.157e-07 * sqMzOverCharge) + 
		        (1.417e-01 * charge));
	}
	
	public static float getCCSFromIMS(float IMS, double precursorMz, byte charge) {
		double gasMass=28.013;
		double temp=31.85;
		double t_diff=273.15;
	    double constant = 18509.8632163405;
	    double reduced_mass = (precursorMz * charge * gasMass) / (precursorMz * charge + gasMass);
	    return (float)((constant * charge) / (Math.sqrt(reduced_mass * (temp + t_diff)) * 1 / IMS));
	}
}
