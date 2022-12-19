package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.ScoredPSM;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public class ScoredPSMFilter implements ScoredPSMFilterInterface {
	private final RetentionTimeFilter rtFilter;
	private final MassErrorFilter precursorFilter;
	private final MassErrorFilter fragmentFilter;

	public ScoredPSMFilter(SearchParameters params, ArrayList<ScoredPSM> passingPSMs) {
		ArrayList<XYPoint> rtPoints=new ArrayList<XYPoint>();
		ArrayList<XYPoint> precursorPoints=new ArrayList<XYPoint>();
		ArrayList<XYPoint> fragmentPoints=new ArrayList<XYPoint>();
		
		for (ScoredPSM psm : passingPSMs) {
			String peptideModSeq=psm.getLibraryEntry().getPeptideModSeq();
			boolean isDecoy=psm.getLibraryEntry().isDecoy();
			float acquiredRT=psm.getMSMS().getScanStartTime()/60f;

			rtPoints.add(psm.getRTData()); 

			if (psm.getDeltaPrecursorMass()!=params.getPrecursorTolerance().getToleranceThreshold()) {
				precursorPoints.add(new PeptideXYPoint(acquiredRT, psm.getDeltaPrecursorMass(), isDecoy, peptideModSeq));
			}
			if (psm.getDeltaPrecursorMass()!=params.getFragmentTolerance().getToleranceThreshold()) {
				fragmentPoints.add(new PeptideXYPoint(acquiredRT, psm.getDeltaFragmentMass(), isDecoy, peptideModSeq));
			}
		}
		
		rtFilter=RetentionTimeFilter.getFilter(rtPoints);
		precursorFilter=MassErrorFilter.getFilter(params.getPrecursorTolerance(), 1, precursorPoints);
		fragmentFilter=MassErrorFilter.getFilter(params.getFragmentTolerance(), 2, fragmentPoints);
	}
	
	@Override
	public float getYRT(float xrt) {
		return rtFilter.getYValue(xrt);
	}
	
	@Override
	public void makePlots(ArrayList<ScoredPSM> psms, Optional<File> saveFileSeed) {
		ArrayList<XYPoint> ms1Errors=new ArrayList<>();
		for (ScoredPSM psm : psms) {
			ms1Errors.add(new XYPoint(psm.getMSMS().getScanStartTime()/60f, psm.getDeltaPrecursorMass()));
		}
		precursorFilter.plot(ms1Errors, saveFileSeed);	
		
		ArrayList<XYPoint> ms2Errors=new ArrayList<>();
		for (ScoredPSM psm : psms) {
			ms2Errors.add(new XYPoint(psm.getMSMS().getScanStartTime()/60f, psm.getDeltaFragmentMass()));
		}
		fragmentFilter.plot(ms2Errors, saveFileSeed);	
		
		ArrayList<XYPoint> rts=new ArrayList<>();
		for (ScoredPSM psm : psms) {
			rts.add(psm.getRTData());
		}
		rtFilter.plot(rts, saveFileSeed);	
	}
	
	@Override
	public boolean passesFilter(ScoredPSM psm) {
		float modelRT=psm.getLibraryEntry().getScanStartTime()/60f;
		float actualRT=psm.getMSMS().getScanStartTime()/60f;
		boolean passes=rtFilter.getProbabilityFitsModel(actualRT, modelRT)>=AbstractRetentionTimeFilter.rejectionPValue;
		
		return passes;
	}
	
	@Override
	public float[] getAdditionalScores(ScoredPSM psm) {
		float modelRT=psm.getLibraryEntry().getScanStartTime()/60f;
		float actualRT=psm.getMSMS().getScanStartTime()/60f;

		float deltaRT=Math.abs(rtFilter.getDelta(actualRT, modelRT));
		float deltaPrecursor=precursorFilter.getCorrectedMassError(actualRT, psm.getDeltaPrecursorMass());
		float deltaFragment=fragmentFilter.getCorrectedMassError(actualRT, psm.getDeltaFragmentMass());
		return new float[] {deltaRT, deltaPrecursor, deltaFragment};
	}

	public RetentionTimeFilter getRtFilter() {
		return rtFilter;
	}

	public MassErrorFilter getPrecursorFilter() {
		return precursorFilter;
	}

	public MassErrorFilter getFragmentFilter() {
		return fragmentFilter;
	}
}
