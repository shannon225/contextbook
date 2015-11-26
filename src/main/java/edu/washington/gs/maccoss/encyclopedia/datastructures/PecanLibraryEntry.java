package edu.washington.gs.maccoss.encyclopedia.datastructures;

//@Immutable
public class PecanLibraryEntry extends LibraryEntry {
	private final float euclidianDistance;
	private final boolean isDecoy;

	public PecanLibraryEntry(double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, boolean isDecoy,
			float euclidianDistance) {
		super(precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray);
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
