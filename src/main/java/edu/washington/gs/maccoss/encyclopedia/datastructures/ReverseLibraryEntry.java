package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

//@Immutable
public class ReverseLibraryEntry extends LibraryEntry {
	public ReverseLibraryEntry(String filename, HashSet<String> accessions, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, float[] correlationArray, boolean[] quantifiedIonsArray, AminoAcidConstants aaConstants) {
		super(filename, accessions, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray, quantifiedIonsArray, aaConstants);
	}
	
	@Override
	public boolean isDecoy() {
		return true;
	}
}
