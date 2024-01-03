package edu.washington.gs.maccoss.encyclopedia.algorithms.curve;

public interface AbstractDilutionCurveFittingParameters {

	public float getWindowInMin();
	
	public float getWindowInMin(float rtInSec);

	public int getNumberOfRTAnchors();

	public int getMaxNumberPeptidesPerProtein();

	public int getTargetTotalNumberOfPeptides();

	public float getMinCVForAnchors();

	public float getMinCVForBadAnchors();

	public int getAssayMaxDensity();

	//public String getTargetAccessionNumberKeyword();
	
	public boolean isTargetedProtein(String accession);
	
	public boolean isEliminatedPeptide(String peptideModSeq);

	public boolean isRequireAlignmentRT();
	
	public boolean isUseLineNoise();
}