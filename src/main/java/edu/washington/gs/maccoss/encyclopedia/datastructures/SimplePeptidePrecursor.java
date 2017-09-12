package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class SimplePeptidePrecursor implements PeptidePrecursor {
	private final byte precursorCharge;
	private final String peptideModSeq;
	private final String massCorrectedPeptideModSeq;

	public SimplePeptidePrecursor(String peptideModSeq, byte precursorCharge) {
		this.peptideModSeq=peptideModSeq;
		this.massCorrectedPeptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq);
		this.precursorCharge=precursorCharge;
	}

	@Override
	public String getMassCorrectedPeptideModSeq() {
		return massCorrectedPeptideModSeq;
	}

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public String getPeptideModSeq() {
		return peptideModSeq;
	}
	
	public String getPeptideSeq() {
		return PeptideUtils.getPeptideSeq(peptideModSeq);
	}

	@Override
	public int compareTo(PeptidePrecursor o) {
		int c=getMassCorrectedPeptideModSeq().compareTo(o.getMassCorrectedPeptideModSeq());
		if (c!=0) return c;
		return precursorCharge-o.getPrecursorCharge();
	}
}
