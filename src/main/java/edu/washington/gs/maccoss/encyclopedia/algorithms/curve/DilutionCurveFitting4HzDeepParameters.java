package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

public class DilutionCurveFitting4HzDeepParameters implements AbstractDilutionCurveFittingParameters {
	final int numberOfRTAnchors=10;
	final int maxNumberPeptidesPerProtein=3;
	final int targetTotalNumberOfPeptides=300; // remember to subtract off anchors (total is 160 peptides)
	final float windowInMin=5f; // in minutes!
	final float minCVForAnchors=0.05f;
	final float minCVForBadAnchors=0.75f;
	final int assayMaxDensity=10;
	final String targetAccessionNumberKeyword="HCMV";
	final boolean requireAlignmentRT=true; // turn off for fitting against PRM

	public float getWindowInMin() {
		return windowInMin;
	}

	public int getNumberOfRTAnchors() {
		return numberOfRTAnchors;
	}

	public int getMaxNumberPeptidesPerProtein() {
		return maxNumberPeptidesPerProtein;
	}

	public int getTargetTotalNumberOfPeptides() {
		return targetTotalNumberOfPeptides;
	}

	public float getMinCVForAnchors() {
		return minCVForAnchors;
	}

	public float getMinCVForBadAnchors() {
		return minCVForBadAnchors;
	}

	public int getAssayMaxDensity() {
		return assayMaxDensity;
	}

	public String getTargetAccessionNumberKeyword() {
		return targetAccessionNumberKeyword;
	}

	public boolean isRequireAlignmentRT() {
		return requireAlignmentRT;
	}
}
