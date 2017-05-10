package edu.washington.gs.maccoss.encyclopedia.datastructures;

import gnu.trove.map.hash.TIntIntHashMap;

public interface FastaEntryInterface extends Comparable<FastaEntryInterface> {

	String getAccession();
	
	String getAnnotation();

	String getFilename();

	String getSequence();
	
	public FastaPeptideEntry getSubEntry(String subSequence);
	
	public FastaPeptideEntry getEntryAsPeptide();

	void addStatistics(TIntIntHashMap map);

}