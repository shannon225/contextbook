package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;

public class IntegratedLibraryEntry extends LibraryEntry {
	private final TransitionRefinementData refinementData;
	
	public IntegratedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, TransitionRefinementData refinementData) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray);
		this.refinementData=refinementData;
	}
	
	public Range getRtRange() {
		return refinementData.getRange();
	}
	
	public TransitionRefinementData getRefinementData() {
		return refinementData;
	}
}
