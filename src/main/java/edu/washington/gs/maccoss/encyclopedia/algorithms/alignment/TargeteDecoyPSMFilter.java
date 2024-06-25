package edu.washington.gs.maccoss.encyclopedia.algorithms.alignment;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.ScoredPSM;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYPoint;

public class TargeteDecoyPSMFilter implements ScoredPSMFilterInterface {
	private final RetentionTimeTargetDecoyFilter rtFilter;
	private final MassErrorFilter precursorFilter;
	private final MassErrorFilter fragmentFilter;
	private final SearchParameters params;

	public TargeteDecoyPSMFilter(SearchParameters params, ArrayList<PSMInterface> passingPSMs, ArrayList<PSMInterface> decoyPSMs) {
		this.params=params;
		
		ArrayList<XYPoint> rtPoints=new ArrayList<XYPoint>();
		ArrayList<XYPoint> precursorPoints=new ArrayList<XYPoint>();
		ArrayList<XYPoint> fragmentPoints=new ArrayList<XYPoint>();
		
		for (PSMInterface psm : passingPSMs) {
			String peptideModSeq=psm.getLibraryEntry().getPeptideModSeq();
			boolean isDecoy=psm.getRTData().isDecoy();
			float acquiredRT=(float)psm.getRTData().y;

			if (!isDecoy) {
				rtPoints.add(psm.getRTData());
			}
			if (params.getPrecursorTolerance().getToleranceThreshold()!=psm.getDeltaPrecursorMass()) {
				precursorPoints.add(new PeptideXYPoint(acquiredRT, psm.getDeltaPrecursorMass(), isDecoy, peptideModSeq));
			}
			if (params.getFragmentTolerance().getToleranceThreshold()!=psm.getDeltaFragmentMass()) {
				fragmentPoints.add(new PeptideXYPoint(acquiredRT, psm.getDeltaFragmentMass(), isDecoy, peptideModSeq));
			}
		}

		ArrayList<XYPoint> decoyRTPoints=new ArrayList<XYPoint>();
		for (PSMInterface psm : decoyPSMs) {
			decoyRTPoints.add(psm.getRTData()); 
		}
		
		rtFilter=RetentionTimeTargetDecoyFilter.getFilter(rtPoints, decoyRTPoints);
		precursorFilter=MassErrorFilter.getFilter(params.getPrecursorTolerance(), 1, precursorPoints);
		fragmentFilter=MassErrorFilter.getFilter(params.getFragmentTolerance(), 2, fragmentPoints);
	}
	
	@Override
	public float getYRT(float xrt) {
		return rtFilter.getYValue(xrt);
	}
	
	@Override
	public void makePlots(SearchParameters params, ArrayList<PSMInterface> psms, Optional<File> saveFileSeed) {
		ArrayList<XYPoint> ms1Errors=new ArrayList<>();
		for (PSMInterface psm : psms) {
			if (params.getPrecursorTolerance().getToleranceThreshold()!=psm.getDeltaPrecursorMass()) {
				ms1Errors.add(new XYPoint(psm.getRTData().y, psm.getDeltaPrecursorMass()));
			}
		}
		precursorFilter.plot(ms1Errors, saveFileSeed);	
		
		ArrayList<XYPoint> ms2Errors=new ArrayList<>();
		for (PSMInterface psm : psms) {
			if (params.getFragmentTolerance().getToleranceThreshold()!=psm.getDeltaFragmentMass()) {
				ms2Errors.add(new XYPoint(psm.getRTData().y, psm.getDeltaFragmentMass()));
			}
		}
		fragmentFilter.plot(ms2Errors, saveFileSeed);	
		
		ArrayList<XYPoint> rts=new ArrayList<>();
		for (PSMInterface psm : psms) {
			rts.add(psm.getRTData());
		}
		rtFilter.plot(rts, saveFileSeed);	
	}
	
	@Override
	public boolean passesFilter(PSMInterface psm) {
		float modelRT=(float)psm.getRTData().x;
		float actualRT=(float)psm.getRTData().y;
		boolean passes=rtFilter.getProbabilityFitsModel(actualRT, modelRT)>=rtFilter.getRejectionPValue();
		
		return passes;
	}
	
	public boolean passesFilter(ScoredPSM psm, float rejectionPValue) {
		float modelRT=psm.getLibraryEntry().getScanStartTime()/60f;
		float actualRT=psm.getMSMS().getScanStartTime()/60f;
		boolean passes=rtFilter.getProbabilityFitsModel(actualRT, modelRT)>=rejectionPValue;
		
		return passes;
	}
	
	@Override
	public float[] getAdditionalScores(PSMInterface psm) {
		float modelRT=(float)psm.getRTData().x;
		float actualRT=(float)psm.getRTData().y;
		
		float deltaRT=Math.abs(rtFilter.getDelta(actualRT, modelRT));
		float deltaPrecursor;
		if (params.getPrecursorTolerance().getToleranceThreshold()!=psm.getDeltaPrecursorMass()) {
			deltaPrecursor=precursorFilter.getCorrectedMassError(actualRT, psm.getDeltaPrecursorMass());
		} else {
			deltaPrecursor=psm.getDeltaPrecursorMass();
		}
		float deltaFragment;
		if (params.getFragmentTolerance().getToleranceThreshold()!=psm.getDeltaFragmentMass()) {
			deltaFragment=fragmentFilter.getCorrectedMassError(actualRT, psm.getDeltaFragmentMass());
		} else {
			deltaFragment=psm.getDeltaFragmentMass();
		}
		
		return new float[] {deltaRT, deltaPrecursor, deltaFragment};
	}

	public RetentionTimeTargetDecoyFilter getRtFilter() {
		return rtFilter;
	}

	public MassErrorFilter getPrecursorFilter() {
		return precursorFilter;
	}

	public MassErrorFilter getFragmentFilter() {
		return fragmentFilter;
	}
}
