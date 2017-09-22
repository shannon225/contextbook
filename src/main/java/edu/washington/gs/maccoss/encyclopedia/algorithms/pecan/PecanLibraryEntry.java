package edu.washington.gs.maccoss.encyclopedia.algorithms.pecan;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FastaPeptideEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;

//@Immutable
public class PecanLibraryEntry extends LibraryEntry {
	private final FastaPeptideEntry entry;
	private final float euclidianDistance;
	private final boolean isDecoy;

	public PecanLibraryEntry(FastaPeptideEntry entry, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, boolean isDecoy,
			float euclidianDistance) {
		super(entry.getFilename(), entry.getAccessions(), precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray);
		this.entry=entry;
		this.isDecoy=isDecoy;
		this.euclidianDistance=euclidianDistance;
	}

	@Override
	public boolean isDecoy() {
		return isDecoy;
	}

	public float getEuclidianDistance() {
		return euclidianDistance;
	}

	public float[] getUnnormalizedIntensities() {
		float[] intensities=getIntensityArray().clone();
		for (int i=0; i<intensities.length; i++) {
			intensities[i]=intensities[i]*euclidianDistance;
		}
		return intensities;
	}
}
