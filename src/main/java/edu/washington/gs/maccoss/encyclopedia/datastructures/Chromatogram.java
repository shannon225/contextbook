package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public interface Chromatogram extends Spectrum {
	public Range getRtRange();
	public float[] getCorrelationArray();
	public float[] getMedianChromatogram();
}
