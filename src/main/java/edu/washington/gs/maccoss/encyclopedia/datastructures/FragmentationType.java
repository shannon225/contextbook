package edu.washington.gs.maccoss.encyclopedia.datastructures;

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
}
