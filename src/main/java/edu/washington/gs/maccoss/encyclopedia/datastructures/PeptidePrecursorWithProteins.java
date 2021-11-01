package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

public interface PeptidePrecursorWithProteins extends PeptidePrecursor {

	public byte getPrecursorCharge();
//	public String getLegacyPeptideModSeq();
	public String getPeptideSeq();
	public String getPeptideModSeq();
	public HashSet<String> getAccessions();
	public float getScore();
}
