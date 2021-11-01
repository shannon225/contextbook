package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

public interface PeptidePrecursor extends Comparable<PeptidePrecursor> {

	public byte getPrecursorCharge();
//	public String getLegacyPeptideModSeq();
	public String getPeptideSeq();
	public String getPeptideModSeq();
}
