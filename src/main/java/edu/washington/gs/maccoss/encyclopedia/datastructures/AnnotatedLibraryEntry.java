package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

public class AnnotatedLibraryEntry extends LibraryEntry {
	private final String[] ionAnnotations;
	
	public AnnotatedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, float[] correlationArray, String[] ionAnnotations) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray);
		this.ionAnnotations=ionAnnotations;
	}
	
	public String[] getIonAnnotations() {
		return ionAnnotations;
	}
}
