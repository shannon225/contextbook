package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;

public interface PeakLocationInferrerInterface {

	Optional<QuantitativeDIAData> getQuantitativeData(TransitionRefinementData data);

	double[] getTopNBestIons(String peptideModSeq, byte precursorCharge);

	float getPreciseRTInSec(SearchJobData job, String peptideModSeq, float detectedRTInSec);

	float getWarpedRTInSec(SearchJobData job, String peptideModSeq);

	List<RetentionTimeAlignmentInterface.AlignmentDataPoint> getAlignmentData(SearchJobData job);
}