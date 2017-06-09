package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;

public interface PeakRTLocatorInterface {

	Optional<Pair<Float, Integer>> getTopNIntensity(TransitionRefinementData data);

	double[] getTopNBestIons(String peptideModSeq);

	float getPreciseRTInSec(SearchJobData job, String peptideModSeq, float detectedRTInSec);

	float getWarpedRTInSec(SearchJobData job, String peptideModSeq);

}