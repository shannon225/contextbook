package edu.washington.gs.maccoss.encyclopedia.algorithms.percolator;

import java.util.ArrayList;

import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;

public class MProphetResult {
	private final ArrayList<PercolatorPeptide> passingPeptides;
	private final LinearDiscriminantAnalysis lda;
	private final ArrayList<String> featureNames;
	private final float pi0;
	
	public MProphetResult(ArrayList<PercolatorPeptide> passingPeptides, LinearDiscriminantAnalysis lda, ArrayList<String> featureNames, float pi0) {
		this.passingPeptides = passingPeptides;
		this.lda = lda;
		this.featureNames=featureNames;
		this.pi0 = pi0;
	}
	
	public ArrayList<PercolatorPeptide> getPassingPeptides() {
		return passingPeptides;
	}
	public LinearDiscriminantAnalysis getLDA() {
		return lda;
	}
	public ArrayList<String> getFeatureNames() {
		return featureNames;
	}
	public float getPi0() {
		return pi0;
	}
	
}
