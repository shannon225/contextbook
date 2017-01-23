package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public enum MassErrorUnitType {
	PPM, AMU;
	public static String toString(MassErrorUnitType type) {
		switch (type) {
			case PPM:
				return "PPM";
			case AMU:
				return "AMU";
			default:
				return "Unknown";
		}
	}
	
	public static MassErrorUnitType getUnitType(String s) {
		if ("PPM".equalsIgnoreCase(s)) return PPM;
		if ("AMU".equalsIgnoreCase(s)) return AMU;
		return null;
	}
}
