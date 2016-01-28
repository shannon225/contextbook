package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPSM;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.Function;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TFloatArrayList;

public class RetentionTimeFilter {
	private final Function rtWarper;
	public RetentionTimeFilter(Function rtWarper) {
		this.rtWarper=rtWarper;
	}
	
	public PercolatorData filterData(PercolatorData perc, ArrayList<PeptideScoringResult> data, SearchParameters parameters) {
		ArrayList<PercolatorPSM> psms=perc.getPsms(); //already reverse sorted

		HashSet<String> passingPSMIDs=new HashSet<String>();
		for (PercolatorPSM psm : psms) {
			if (psm.getQValue()<=parameters.getPercolatorThreshold()) {
				passingPSMIDs.add(psm.getPsmID());
			}
		}
		
		TFloatArrayList deltas=new TFloatArrayList();

		for (PeptideScoringResult result : data) {
			if (result.getGoodStripes().size()>0) {
				String peptideModSeq=result.getEntry().getPeptideModSeq();
				if (passingPSMIDs.contains(peptideModSeq+"+"+result.getEntry().getPrecursorCharge())) {
					LibraryEntry entry=result.getEntry();
					float entryTime=rtWarper.getYValue(entry.getRetentionTime());
					
					Pair<ScoredObject<Stripe>, float[]> first=result.getGoodStripes().get(0);
					float deltaRT=first.x.y.getScanStartTime()/60f-entryTime;
					deltas.add(deltaRT);					
				}
			}
		}
		
		float[] deltaArray=deltas.toArray();
		float mean=General.mean(deltaArray);
		float stdev=General.stdev(deltaArray);
		float upperThreshold=mean+2.0f*stdev;
		float lowerThreshold=mean-2.0f*stdev;

		HashSet<String> rtFilteredPSMIDs=new HashSet<String>();
		for (PeptideScoringResult result : data) {
			if (result.getGoodStripes().size()>0) {
				String peptideModSeq=result.getEntry().getPeptideModSeq();
				LibraryEntry entry=result.getEntry();
				float entryTime=rtWarper.getYValue(entry.getRetentionTime());

				Pair<ScoredObject<Stripe>, float[]> first=result.getGoodStripes().get(0);
				float deltaRT=first.x.y.getScanStartTime()/60f-entryTime;

				if (deltaRT<=upperThreshold&&deltaRT>=lowerThreshold) {
					rtFilteredPSMIDs.add(peptideModSeq);
				}
			}
		}

		int decoys=0;
		int nondecoys=0;
		TFloatArrayList qvalues=new TFloatArrayList();
		for (PercolatorPSM psm : psms) {
			if (rtFilteredPSMIDs.contains(psm.getPsmID())) {
				if (psm.isDecoy()) {
					decoys++;
				} else {
					nondecoys++;
				}
			}
			qvalues.add(decoys/(float)(decoys+nondecoys));
		}
		
		// convert FDRs to q-values
		float minValue=Float.MAX_VALUE;
		for (int i=qvalues.size()-1; i>=0; i--) {
			if (qvalues.get(i)>minValue) {
				qvalues.set(i, minValue);
			} else {
				minValue=qvalues.get(i);
			}
		}

		ArrayList<PercolatorPSM> rtFilteredPSMs=new ArrayList<PercolatorPSM>();
		for (int i=0; i<psms.size(); i++) {
			if (rtFilteredPSMIDs.contains(psms.get(i).getPsmID())) {
				rtFilteredPSMs.add(psms.get(i).clone(qvalues.get(i)));
			}
		}
		
		return perc.clone(rtFilteredPSMs);
	}
}
