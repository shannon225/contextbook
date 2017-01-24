package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public class XCorrCalculator {
	public static Spectrum normalize(Spectrum s, double precursorMz) {
		double[] masses=s.getMassArray();
		float[] intensities=s.getIntensityArray();
		
		double minimumPrecursorRemoved=precursorMz-5.0;
		double maximumPrecursorRemoved=precursorMz+5.0;
		
		return s; // FIXME FINISH ALGORITHM!
	}

}
