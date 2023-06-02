package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

public class Targeted10HzParametersWithWideWindows implements AbstractDilutionCurveFittingParameters {
	final int numberOfRTAnchors=0;
	final int maxNumberPeptidesPerProtein=3;
	final int targetTotalNumberOfPeptides=3000; // remember to subtract off anchors (total is 160 peptides)
	final float windowInMin=12f; // in minutes!
	final float minCVForAnchors=0.05f;
	final float minCVForBadAnchors=0.75f;
	final int assayMaxDensity=25;
	final String targetAccessionNumberKeyword="";
	final boolean requireAlignmentRT=true; // turn off for fitting against PRM
	final boolean useLineNoise=false; // newer versions should set this to "true"

	public float getWindowInMin() {
		return windowInMin;
	}
	public float getWindowInMin(float rtInSec) {
		if (rtInSec/60>70) return windowInMin*2f;
		if (rtInSec/60>30) return windowInMin*1.5f;
		if (rtInSec/60<30) return windowInMin*1.25f;
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
	
	@Override
	public boolean isTargetedProtein(String accession) {
		return true;
	}

	public boolean isRequireAlignmentRT() {
		return requireAlignmentRT;
	}

	public boolean isUseLineNoise() {
		return useLineNoise;
	}
}
