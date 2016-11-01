package edu.washington.gs.maccoss.encyclopedia.datastructures;

public class SimplePeptidePrecursor implements PeptidePrecursor {
	private final byte precursorCharge;
	private final String peptideModSeq;

	public SimplePeptidePrecursor(String peptideModSeq, byte precursorCharge) {
		this.peptideModSeq=peptideModSeq;
		this.precursorCharge=precursorCharge;
	}

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public String getPeptideModSeq() {
		return peptideModSeq;
	}

	@Override
	public int compareTo(PeptidePrecursor o) {
		int c=peptideModSeq.compareTo(o.getPeptideModSeq());
		if (c!=0) return c;
		return precursorCharge-o.getPrecursorCharge();
	}
}
