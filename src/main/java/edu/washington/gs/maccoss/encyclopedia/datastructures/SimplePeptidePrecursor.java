package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;

public class SimplePeptidePrecursor implements PeptidePrecursor {
	private final byte precursorCharge;
	private final String peptideModSeq;
	private final String massCorrectedPeptideModSeq;

	public SimplePeptidePrecursor(String peptideModSeq, byte precursorCharge, AminoAcidConstants aaConstants) {
		this.peptideModSeq=peptideModSeq;
		this.massCorrectedPeptideModSeq=PeptideUtils.getCorrectedMasses(peptideModSeq, aaConstants);
		this.precursorCharge=precursorCharge;
	}

	@Override
	public String getPeptideModSeq() {
		return massCorrectedPeptideModSeq;
	}

	public byte getPrecursorCharge() {
		return precursorCharge;
	}

	public String getLegacyPeptideModSeq() {
		return peptideModSeq;
	}
	
	public String getPeptideSeq() {
		return PeptideUtils.getPeptideSeq(peptideModSeq);
	}

	@Override
	public int compareTo(PeptidePrecursor o) {
		int c=getPeptideModSeq().compareTo(o.getPeptideModSeq());
		if (c!=0) return c;
		return precursorCharge-o.getPrecursorCharge();
	}
	
	@Override
	public int hashCode() {
		return getPeptideModSeq().hashCode()+precursorCharge;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof PeptidePrecursor) {
			return compareTo((PeptidePrecursor)obj)==0;
		}
		return false;
	}
}
