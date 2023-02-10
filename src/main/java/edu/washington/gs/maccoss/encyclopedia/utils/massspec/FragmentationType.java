package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public enum FragmentationType {
	CID, ETD, HCD, SILAC;
	public static String toString(FragmentationType type) {
		switch (type) {
			case CID:
				return "CID";
			case ETD:
				return "ETD";
			case HCD:
				return "HCD";
			case SILAC:
				return "SILAC";
			default:
				return "Unknown";
		}
	}
	public static String toName(FragmentationType type) {
		switch (type) {
			case CID:
				return "CID/HCD (B/Y)";
			case ETD:
				return "ETD (C/Z/Z+1)";
			case HCD:
				return "HCD (Y-Only)";
			case SILAC:
				return "CID/HCD SILAC (Y-Only)";
			default:
				return "Unknown";
		}
	}
	
	public static FragmentationType getFragmentationType(String s) {
		if ("CID".equalsIgnoreCase(s)) return CID;
		if ("CID (B/Y)".equalsIgnoreCase(s)) return CID;
		if ("CID/HCD (B/Y)".equalsIgnoreCase(s)) return CID;
		
		if ("HCD".equalsIgnoreCase(s)) return HCD;
		if ("YONLY".equalsIgnoreCase(s)) return HCD;
		if ("HCD (Y-Only)".equalsIgnoreCase(s)) return HCD;
		
		if ("SILAC".equalsIgnoreCase(s)) return SILAC;
		if ("CID/HCD SILAC (Y-Only)".equalsIgnoreCase(s)) return SILAC;
		
		if ("ETD".equalsIgnoreCase(s)) return ETD;
		if ("ETD (C/Z/Z+1)".equalsIgnoreCase(s)) return ETD;
		return null;
	}
}
