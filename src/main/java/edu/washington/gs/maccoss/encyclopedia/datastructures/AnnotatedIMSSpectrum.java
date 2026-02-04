package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import gnu.trove.map.hash.TDoubleObjectHashMap;

public class AnnotatedIMSSpectrum extends AnnotatedSpectrum {

	public AnnotatedIMSSpectrum(Spectrum s, PeptidePrecursor entry, SearchParameters parameters) {
		super(s, entry, parameters);
	}

	public AnnotatedIMSSpectrum(String peptideModSeq, Spectrum s, TDoubleObjectHashMap<String> annotationMap, byte precursorCharge, MassTolerance tolerance) {
		super(peptideModSeq, s, annotationMap, precursorCharge, tolerance);
	}

	@Override
	public GraphType getType() {
		if (getIonMobilityArray().isPresent()) return GraphType.imsspectrum;
		return GraphType.spectrum;
	}
}
