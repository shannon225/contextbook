package edu.washington.gs.maccoss.encyclopedia.datastructures;

public interface PeptidePrecursor extends Comparable<PeptidePrecursor> {

	public byte getPrecursorCharge();
	public String getPeptideModSeq();
	public String getPeptideSeq();
	public String getMassCorrectedPeptideModSeq();

}
