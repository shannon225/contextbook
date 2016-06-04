package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public enum FragmentationType {
	CID, ETD, YONLY;
	public static String toString(FragmentationType type) {
		switch (type) {
			case CID:
				return "CID";
			case ETD:
				return "ETD";
			case YONLY:
				return "YONLY";
			default:
				return "Unknown";
		}
	}
	
	public static FragmentationType getFragmentationType(String s) {
		if ("CID".equalsIgnoreCase(s)) return CID;
		if ("ETD".equalsIgnoreCase(s)) return ETD;
		if ("YONLY".equalsIgnoreCase(s)) return YONLY;
		if ("CID (B/Y)".equalsIgnoreCase(s)) return CID;
		if ("ETD (C/Z/Z+1)".equalsIgnoreCase(s)) return ETD;
		if ("HCD (Y-Only)".equalsIgnoreCase(s)) return YONLY;
		return null;
	}
}
