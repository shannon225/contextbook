package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;

public enum ScoringBreadthType {
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
	
	public static ScoringBreadthType getType(String s) {
		for (ScoringBreadthType type : ScoringBreadthType.values()) {
			if (type.toShortname().equalsIgnoreCase(s)) {
				return type;
			}
		}
		throw new EncyclopediaException("Unexpected scoring breadth type ["+s+"]");
	}
	
	public boolean runRecalibration() {
		switch (this) {
		case ENTIRE_RT_WINDOW:
			return true;
		case RECALIBRATED_20_PERCENT:
			return true;
		case RECALIBRATED_PEAK_WIDTH:
			return true;
		case UNCALIBRATED_20_PERCENT:
			return false;
		case UNCALIBRATED_PEAK_WIDTH:
			return false;
		}
		return true;
	}

	public String toShortname() {
		switch (this) {
		case ENTIRE_RT_WINDOW:
			return "window";
		case RECALIBRATED_20_PERCENT:
			return "recal20";
		case RECALIBRATED_PEAK_WIDTH:
			return "recal";
		case UNCALIBRATED_20_PERCENT:
			return "uncal20";
		case UNCALIBRATED_PEAK_WIDTH:
			return "uncal";
		}
		return "Unknown";
	};
}
