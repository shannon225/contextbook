package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

public interface AbstractDilutionCurveFittingParameters {

	public float getWindowInMin();

	public int getNumberOfRTAnchors();

	public int getMaxNumberPeptidesPerProtein();

	public int getTargetTotalNumberOfPeptides();

	public float getMinCVForAnchors();

	public float getMinCVForBadAnchors();

	public int getAssayMaxDensity();

	public String getTargetAccessionNumberKeyword();

	public boolean isRequireAlignmentRT();
}