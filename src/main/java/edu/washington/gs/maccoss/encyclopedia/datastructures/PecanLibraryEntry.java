package edu.washington.gs.maccoss.encyclopedia.datastructures;

import jdk.nashorn.internal.ir.annotations.Immutable;

@Immutable
public class PecanLibraryEntry extends LibraryEntry {
	public PecanLibraryEntry(double precursorMZ, byte precursorCharge,
			String peptideModSeq, int copies, float retentionTime, float score,
			double[] massArray, float[] intensityArray) {
		super(precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime,
				score, massArray, intensityArray);
	}
}
