package edu.washington.gs.maccoss.encyclopedia.datastructures;

import jdk.nashorn.internal.ir.annotations.Immutable;

@Immutable
public class ReverseLibraryEntry extends LibraryEntry {
	public ReverseLibraryEntry(double precursorMZ, byte precursorCharge,
			String peptideModSeq, int copies, float retentionTime, float score,
			double[] massArray, float[] intensityArray) {
		super(precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime,
				score, massArray, intensityArray);
	}
}
