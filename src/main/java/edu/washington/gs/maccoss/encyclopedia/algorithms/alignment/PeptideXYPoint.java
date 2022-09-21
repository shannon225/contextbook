package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public class PeptideXYPoint extends XYPoint {
	private final boolean isDecoy;
	private final String peptideModSeq;

	public PeptideXYPoint(double x, double y, boolean isDecoy, String peptideModSeq) {
		super(x, y);
		this.isDecoy=isDecoy;
		this.peptideModSeq=peptideModSeq;
	}
	
	public String getPeptideModSeq() {
		return peptideModSeq;
	}
	
	public boolean isDecoy() {
		return isDecoy;
	}
}