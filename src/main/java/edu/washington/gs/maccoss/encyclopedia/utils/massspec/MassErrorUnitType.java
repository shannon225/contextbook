package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public enum MassErrorUnitType {
	PPM, AMU, RESOLUTION;
	public static String toString(MassErrorUnitType type) {
		switch (type) {
			case RESOLUTION:
				return "Resolution";
			case PPM:
				return "PPM";
			case AMU:
				return "AMU";
			default:
				return "Unknown";
		}
	}
	
	public static MassErrorUnitType getUnitType(String s) {
		if ("Resolution".equalsIgnoreCase(s)) return RESOLUTION;
		if ("PPM".equalsIgnoreCase(s)) return PPM;
		if ("AMU".equalsIgnoreCase(s)) return AMU;
		return null;
	}
}
