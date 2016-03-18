package edu.washington.gs.maccoss.encyclopedia.datastructures;

public enum DataAcquisitionType {
	DIA, OVERLAPPING_DIA, DDA;
	public static String toString(DataAcquisitionType type) {
		switch (type) {
			case DDA:
				return "DDA";
			case DIA:
				return "DIA";
			case OVERLAPPING_DIA:
				return "OverlappingDIA";
			default:
				return "Unknown";
		}
	}
	
	public static DataAcquisitionType getAcquisitionType(String s) {
		if ("DDA".equalsIgnoreCase(s)) return DDA;
		if ("DIA".equalsIgnoreCase(s)) return DIA;
		if ("OVERLAPPINGDIA".equalsIgnoreCase(s)) return OVERLAPPING_DIA;
		if ("OVERLAPPING".equalsIgnoreCase(s)) return OVERLAPPING_DIA;
		if ("OVERLAPPING DIA".equalsIgnoreCase(s)) return OVERLAPPING_DIA;
		if ("NON-OVERLAPPING DIA".equalsIgnoreCase(s)) return DIA;
		return null;
	}
}
