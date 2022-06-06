package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

public class Targeted10HzParameters implements AbstractDilutionCurveFittingParameters {
	final int numberOfRTAnchors=0;
	final int maxNumberPeptidesPerProtein=5;
	final int targetTotalNumberOfPeptides=3000; // remember to subtract off anchors (total is 160 peptides)
	final float windowInMin=5f; // in minutes!
	final float minCVForAnchors=0.05f;
	final float minCVForBadAnchors=0.75f;
	final int assayMaxDensity=28;
	final String targetAccessionNumberKeyword="";
	final boolean requireAlignmentRT=true; // turn off for fitting against PRM
	final boolean useLineNoise=false; // newer versions should set this to "true"

	public float getWindowInMin() {
		return windowInMin;
	}
	public float getWindowInMin(float rtInSec) {
		if (rtInSec/60>60) return windowInMin*3f;
		if (rtInSec/60>50) return windowInMin*2f;
		if (rtInSec/60>45) return windowInMin*1.5f;
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

	public boolean isUseLineNoise() {
		return useLineNoise;
	}
}
