package edu.washington.gs.maccoss.encyclopedia.datastructures;

public enum DataAcquisitionType {
	DIA, OVERLAPPING_DIA;
	public static String toString(DataAcquisitionType type) {
		switch (type) {
			case DIA:
				return "DIA";
			case OVERLAPPING_DIA:
				return "OverlappingDIA";
			default:
				return "Unknown";
		}
	}
	
	public static String toName(DataAcquisitionType type) {
		switch (type) {
			case DIA:
				return "Non-Overlapping DIA";
			case OVERLAPPING_DIA:
				return "Overlapping DIA";
			default:
				return "Unknown";
		}
	}
	
	public static DataAcquisitionType getAcquisitionType(String s) {
		if ("FALSE".equalsIgnoreCase(s)) return DIA;
		if ("DIA".equalsIgnoreCase(s)) return DIA;
		if ("TRUE".equalsIgnoreCase(s)) return OVERLAPPING_DIA;
		if ("OVERLAPPINGDIA".equalsIgnoreCase(s)) return OVERLAPPING_DIA;
		if ("OVERLAPPING".equalsIgnoreCase(s)) return OVERLAPPING_DIA;
		if ("OVERLAPPING DIA".equalsIgnoreCase(s)) return OVERLAPPING_DIA;
		if ("NON-OVERLAPPING DIA".equalsIgnoreCase(s)) return DIA;
		return null;
	}
}
