package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.RetentionTimeAlignmentInterface.AlignmentDataPoint;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefinementData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.Quadruplet;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Peak;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakChromatogram;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.QuantitativeDIAData;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.procedure.TObjectFloatProcedure;

public class LibraryPeakLocationInferrer implements PeakLocationInferrerInterface {
	private final SearchParameters params;
	private final HashMap<String, LibraryEntry> entriesByPeptideKey;
	private final HashMap<SearchJobData, TObjectFloatHashMap<String>> rtByPeptideModSeq;
	
	public LibraryPeakLocationInferrer(Collection<LibraryEntry> entries, HashMap<SearchJobData, TObjectFloatHashMap<String>> rtByPeptideModSeq, SearchParameters params) {
		entriesByPeptideKey=new HashMap<String, LibraryEntry>();
		for (LibraryEntry entry : entries) {
			entriesByPeptideKey.put(getPeptideKey(entry.getPeptideModSeq(), entry.getPrecursorCharge()), entry);
		}
		this.params = params;
		this.rtByPeptideModSeq=rtByPeptideModSeq;
	}
	
	public static String getPeptideKey(String peptideModSeq, byte precursorCharge) {
		return peptideModSeq+"_"+precursorCharge;
	}

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
			return Optional.of(new QuantitativeDIAData(data.getPeptideModSeq(), data.getPrecursorCharge(), data.getApexRT(), data.getRange(), topNMasses, topNIntensities, topNCorrelations, data.getIonMobility(), params.getAAConstants()));
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
		return Optional.of(new QuantitativeDIAData(data.getPeptideModSeq(), data.getPrecursorCharge(), data.getApexRT(), data.getRange(), topNMasses, topNIntensities, topNCorrelations, data.getIonMobility(), params.getAAConstants()));

	}

	@Override
	public double[] getTopNBestIons(String peptideModSeq, byte precursorCharge) {
		LibraryEntry entry=entriesByPeptideKey.get(getPeptideKey(peptideModSeq, precursorCharge));
		ArrayList<Peak> peaks=entry.getPeaksByIntensityFraction(0.1f, params.getNumberOfQuantitativePeaks());
		
		double[] masses=new double[peaks.size()];
		for (int i = 0; i < masses.length; i++) {
			masses[i]=peaks.get(i).mass;
		}
		
		return masses;
	}

	@Override
	public float getPreciseRTInSec(SearchJobData job, String peptideModSeq, float detectedRTInSec) {
		return rtByPeptideModSeq.get(job).get(peptideModSeq);
	}

	@Override
	public float getWarpedRTInSec(SearchJobData job, String peptideModSeq) {
		return rtByPeptideModSeq.get(job).get(peptideModSeq);
	}

	@Override
	public List<AlignmentDataPoint> getAlignmentData(SearchJobData job) {
		TObjectFloatHashMap<String> map=rtByPeptideModSeq.get(job);
		
		ArrayList<AlignmentDataPoint> data=new ArrayList<RetentionTimeAlignmentInterface.AlignmentDataPoint>();
		map.forEachEntry(new TObjectFloatProcedure<String>() {
			@Override
			public boolean execute(String a, float b) {
				LibraryEntry entry=entriesByPeptideKey.get(a);
				if (entry!=null) {
					data.add(AlignmentDataPoint.of(entry.getRetentionTimeInSec()/60f, b/60f, b/60f, 0, 0, false, a));
				}
				return true;
			}
		});
		
		return data;
	}
	

}
