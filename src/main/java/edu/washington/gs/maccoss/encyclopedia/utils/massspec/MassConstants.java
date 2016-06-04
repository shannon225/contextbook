package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public class MassConstants {
	public final static double neutronMass=1.0086649158849;
	public final static double protonMass=1.00727646681290;
	public final static double hydrogenMass=1.007825032071;
	public final static double oh2=15.9949146195616+2*hydrogenMass;
	
	public static Float getModificationMass(String mod) {
		if ("Cam".equals(mod)) {
			return 57.0214635f;
		} else if ("O".equals(mod)) {
			return 15.994915f;
		}
		return null;
	}
	
	private static final MassTolerance tolerance=new MassTolerance(1.0); // 1 ppm is about the accuracy of floats 
	public static double getNeutralLoss(double modificationMass) {
		if (isPhosphoMass(modificationMass)) {
			return 97.976896;
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
