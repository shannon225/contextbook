package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public class RTRTPoint extends XYPoint {
	private final boolean isDecoy;

	public RTRTPoint(double x, double y, boolean isDecoy) {
		super(x, y);
		this.isDecoy=isDecoy;
	}
	
	public boolean isDecoy() {
		return isDecoy;
	}
}