package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Quadruplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakChromatogram;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;

public class SimplePeakLocationInferrer implements PeakLocationInferrerInterface {
	private static float RT_OUTLIER_REJECTION_PROBABILITY=0.001f;
	
	// alignments are seed (x) to sample (y), in minutes
	private final HashMap<SearchJobData, RetentionTimeAlignmentInterface> alignmentMap;
	private final HashMap<SearchJobData, List<RetentionTimeAlignmentInterface.AlignmentDataPoint>> alignmentDataMap;

	// alignedRTs are as if they were in the seed (x) file
	private final HashMap<String, Float> alignedRTInMinBySequenceMap;
	
	private final HashMap<String, double[]> bestIons;
	private final SearchParameters params;

	public SimplePeakLocationInferrer(HashMap<SearchJobData, RetentionTimeAlignmentInterface> alignmentMap, HashMap<SearchJobData, List<RetentionTimeAlignmentInterface.AlignmentDataPoint>> alignmentDataMap, HashMap<String, Float> alignedRTInMinBySequenceMap, HashMap<String, double[]> bestIons, SearchParameters params) {
		this.alignmentMap=alignmentMap;
		this.alignmentDataMap = alignmentDataMap;
		this.alignedRTInMinBySequenceMap=alignedRTInMinBySequenceMap;
		this.bestIons=bestIons;
		this.params=params;
	}
	
	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface#getQuantitativeData(edu.washington.gs.maccoss.encyclopedia.algorithms.TransitionRefinementData)
	 */
	@Override
	public Optional<QuantitativeDIAData> getQuantitativeData(TransitionRefinementData data) {
		String peptideModSeq=data.getPeptideModSeq();
		double[] topNMasses=getTopNBestIons(peptideModSeq, data.getPrecursorCharge());
		double[] masses=FragmentIon.getMasses(data.getFragmentMassArray());
		float[] intensities=data.getIntegrationArray();
		float[] correlations=data.getCorrelationArray();

		if (params.getMinNumOfQuantitativePeaks()>0) {
			if (topNMasses==null||topNMasses.length<params.getMinNumOfQuantitativePeaks()) {
				return Optional.empty();
			}
		}
		
		if (topNMasses==null||topNMasses.length==0) {
			ArrayList<PeakChromatogram> topN=data.getTopNPeaks(TransitionRefiner.quantitativeCorrelationThreshold, params.getEffectiveNumberOfQuantitativePeaks());
			Quadruplet<double[], float[], float[], boolean[]> pair=PeakChromatogram.toChromatogramArrays(topN);
			topNMasses=pair.x;
			float[] topNIntensities=pair.y;
			float[] topNCorrelations=pair.z;
			return Optional.of(new QuantitativeDIAData(data.getPeptideModSeq(), data.getPrecursorCharge(), data.getApexRT(), data.getRange(), topNMasses, topNIntensities, topNCorrelations, params.getAAConstants()));
		}
		
		float[] topNIntensities=new float[topNMasses.length];
		float[] topNCorrelations=new float[topNMasses.length];
		for (int i=0; i<topNMasses.length; i++) {
			float sum=0.0f;
			float maxCorr=0.0f;
			int[] optionalIndex=params.getFragmentTolerance().getIndicies(masses, topNMasses[i]);
			for (int index : optionalIndex) {
				sum+=intensities[index];
				maxCorr=Math.max(maxCorr, correlations[index]);
			}
			topNIntensities[i]=sum;
			topNCorrelations[i]=maxCorr;
		}
		return Optional.of(new QuantitativeDIAData(data.getPeptideModSeq(), data.getPrecursorCharge(), data.getApexRT(), data.getRange(), topNMasses, topNIntensities, topNCorrelations, params.getAAConstants()));
	}

	/* (non-Javadoc)
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface#getTopNBestIons(java.lang.String, byte)
	 */
	@Override
	public double[] getTopNBestIons(String peptideModSeq, byte precursorCharge) {
		return bestIons.get(peptideModSeq);
	}
	
	/**
	 * Prefers the detectedRTInSec if it's available (and within the probability model tolerance), otherwise, uses a warped RT
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface#getPreciseRTInSec(edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData, java.lang.String, float)
	 */
	@Override
	public float getPreciseRTInSec(SearchJobData job, String peptideModSeq, float detectedRTInSec) {
		RetentionTimeAlignmentInterface f=alignmentMap.get(job);
		Float alignedRTInMin=alignedRTInMinBySequenceMap.get(peptideModSeq);
		if (alignedRTInMin==null) {
			return detectedRTInSec;
		}
		
		if (f==null) {
			return detectedRTInSec;
		} else {
			float warpedRTInMin=f.getYValue(alignedRTInMin);
			final float actualRT=detectedRTInSec/60f;
			float prob=f.getProbabilityFitsModel(actualRT, actualRT-warpedRTInMin);
			if (prob>RT_OUTLIER_REJECTION_PROBABILITY) {
				return detectedRTInSec;
			} else {
				return warpedRTInMin*60f;
			}
		}
	}

	/**
	 * Always reports the warped RT
	 * @see edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface#getWarpedRTInSec(edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData, java.lang.String)
	 */
	@Override
	public float getWarpedRTInSec(SearchJobData job, String peptideModSeq) {
		RetentionTimeAlignmentInterface f=alignmentMap.get(job);
		Float alignedRTInMin=alignedRTInMinBySequenceMap.get(peptideModSeq);
		if (alignedRTInMin==null) {
			Logger.errorLine("Couldn't find retention time for peptide ("+peptideModSeq+") in file ("+job.getOriginalDiaFileName()+").");
			return -1;
		}
		
		if (f==null) {
			// job is the seed
			return alignedRTInMin*60f;
		} else {
			return f.getYValue(alignedRTInMin)*60f;
		}
	}

	@Override
	public List<RetentionTimeAlignmentInterface.AlignmentDataPoint> getAlignmentData(SearchJobData job) {
		return alignmentDataMap.getOrDefault(job, Collections.emptyList());
	}

	public RetentionTimeAlignmentInterface getAlignment(SearchJobData jobData) {
		return alignmentMap.get(jobData);
	}
}
