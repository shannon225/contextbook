package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public class AnnotatedIMSSpectrum extends AnnotatedSpectrum {

	public AnnotatedIMSSpectrum(Spectrum s, PeptidePrecursor entry, SearchParameters parameters) {
		super(s, entry, parameters);
	}

	public AnnotatedIMSSpectrum(Spectrum s, TDoubleObjectHashMap<String> annotationMap, MassTolerance tolerance) {
		super(s, annotationMap, tolerance);
	}

	public AnnotatedIMSSpectrum(String name, double mz, float scanStartTime, double[] masses, float[] intensities, Optional<float[]> ionMobilityArray,
			FragmentIon[] annotations) {
		super(name, mz, scanStartTime, masses, intensities, ionMobilityArray, annotations);
	}

	@Override
	public GraphType getType() {
		if (getIonMobilityArray().isPresent()) return GraphType.imsspectrum;
		return GraphType.spectrum;
	}
}
