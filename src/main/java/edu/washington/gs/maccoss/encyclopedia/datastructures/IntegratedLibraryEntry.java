package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;

public class IntegratedLibraryEntry extends ChromatogramLibraryEntry {
	private final TransitionRefinementData refinementData;
	
	public IntegratedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, TransitionRefinementData refinementData) {
		this(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, refinementData.getCorrelationArray(), refinementData);
	}

	public IntegratedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, float[] correlationArray, TransitionRefinementData refinementData) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray, refinementData.getMedianChromatogram(), refinementData.getRange(), refinementData.getAaConstants());
		this.refinementData=refinementData;
	}
	
	public TransitionRefinementData getRefinementData() {
		return refinementData;
	}
}
