package edu.washington.gs.maccoss.encyclopedia.algorithms.precursor;

import java.util.ArrayList;
import java.util.Collections;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMScoredSpectrum;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PeakTrace;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.RetentionTimeComparator;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakWithTime;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PrecursorIon;

// mutable! Not thread safe!
public class PrecursorIonMap {
	private final float rtHalfWindowSize;
	private final MassTolerance tolerance;
	private final ArrayList<IntegratedPeptide> peptides=new ArrayList<>();
	int currentIndex=0;
	
	public PrecursorIonMap(MassTolerance tolerance, float rtTotalWindowSize) {
		this.tolerance = tolerance;
		this.rtHalfWindowSize=rtTotalWindowSize/2.0f;
	}

	public void addPeptide(String peptideModSeq, String sourceFile, ArrayList<PSMScoredSpectrum> psms) {
		peptides.add(new IntegratedPeptide(peptideModSeq, sourceFile, psms));
	}
	
	public void finalizeMap() {
		Collections.sort(peptides, RetentionTimeComparator.comparator);
	}
	
	public void processSpectrum(PrecursorScan scan) {
		if (currentIndex>=peptides.size()) {
			return;
		}
		
		float minRT=scan.getScanStartTime()-rtHalfWindowSize;
		float maxRT=scan.getScanStartTime()+rtHalfWindowSize;
		
		while (true) {
			if (currentIndex>=peptides.size()) {
				// finished all peptides, so quit
				return;
			}
			if (peptides.get(currentIndex).getRetentionTimeInSec()>maxRT) {
				// the leading peptide is later than the scanning window, so quit
				return;
			}
			
			if (peptides.get(currentIndex).getRetentionTimeInSec()<minRT) {
				// the leading peptide is before the scanning window, so stop considering it and look at the next one
				currentIndex++;
				continue;
			}
			
			for (int i = currentIndex; i < peptides.size(); i++) {
				IntegratedPeptide pep=peptides.get(i);
				if (pep.getRetentionTimeInSec()>maxRT) {
					// beyond the current window, so quit
					return;
				}
				PrecursorIon[] ions=pep.getIons();
				double[] mzs=new double[ions.length];
				for (int j = 0; j < mzs.length; j++) {
					mzs[j]=ions[j].getMass();
				}
				// all peptides trapped here are in the scanning window
				Pair<double[], float[]> integrations=tolerance.getIntegratedIntensitiesAndMasses(scan.getMassArray(), scan.getIntensityArray(), mzs);
				double[] masses=integrations.x;
				float[] intensities=integrations.y;
				
				for (int j = 0; j < ions.length; j++) {
					pep.getPeaks()[j].add(new PeakWithTime(masses[j], intensities[j], scan.getScanStartTime()));
				}
			}
			
			// processed all peptides against this scan, so quit
			return;
		}
	}
	
	public ArrayList<PeakTrace<PrecursorIon>> getCenteredTraces() {
		ArrayList<PeakTrace<PrecursorIon>> traces=new ArrayList<>();
		for (IntegratedPeptide peptide : peptides) {
			ArrayList<PeakTrace<PrecursorIon>> centeredTraces = peptide.getCenteredTraces();
			for (PeakTrace<PrecursorIon> trace : centeredTraces) {
				if (trace.getIon().getIsotope()>=0) {
					traces.add(trace);
				}
			}
		}
		return traces;
	}
	
	public ArrayList<IntegratedPeptide> getPeptides() {
		return peptides;
	}
}
