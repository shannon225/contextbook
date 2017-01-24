package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

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
	
	public static Float getModificationMass(String mod) {
		if ("Cam".equals(mod)) {
			return 57.0214635f;
		} else if ("O".equals(mod)) {
			return 15.994915f;
		}
		return null;
	}
	
	private static final MassTolerance tolerance=new MassTolerance(1.0); // 1 ppm is about the accuracy of floats 
	public static double getNeutralLoss(char aa, double modificationMass) {
		if (isPhosphoMass(modificationMass)) {
			if (aa=='S'||aa=='T') {
				return 97.976896;
			}
		}
		return 0.0;
	}
	public static boolean isPhosphoMass(double modificationMass) {
		if (tolerance.equals(80.0, modificationMass)) {
			return true;
		} else if (tolerance.equals(79.966331, modificationMass)) {
			return true;
		}
		return false;
	}
}
