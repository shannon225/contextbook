package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

public enum CAPSiLScoringBreadthType {
	ENTIRE_RT_WINDOW, RECALIBRATED_20_PERCENT, RECALIBRATED_PEAK_WIDTH, UNCALIBRATED_20_PERCENT, UNCALIBRATED_PEAK_WIDTH;
	
	public String toString() {
		switch (this) {
		case ENTIRE_RT_WINDOW:
			return "Across entire window";
		case RECALIBRATED_20_PERCENT:
			return "Recalibrated (20% gradient)";
		case RECALIBRATED_PEAK_WIDTH:
			return "Recalibrated (peak width only)";
		case UNCALIBRATED_20_PERCENT:
			return "Uncalibrated (20% gradient)";
		case UNCALIBRATED_PEAK_WIDTH:
			return "Uncalibrated (peak width only)";
		}
		return "Unknown";
	};
}
