package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

//@Immutable
public class SearchParameters {
	private final AminoAcidConstants aaConstants;
	private final FragmentationType fragType;
	private final MassTolerance precursorTolerance;
	private final MassTolerance fragmentTolerance;
	private final DigestionEnzyme enzyme;
	private final int minPeptideLength;
	private final int maxPeptideLength;
	private final int maxMissedCleavages;
	private final byte minCharge;
	private final byte maxCharge;
	private final int minEluteTime;
	private final int numberOfReportedPeaks;
	private final boolean addDecoysToBackgound;
	private final boolean dontRunDecoys; // only for testing
	private final float percolatorThreshold;
	
	public String toString() {
		final StringBuilder sb=new StringBuilder();
		sb.append(" -fixed "+aaConstants.getFixedModString()+"\n");
		sb.append(" -frag "+FragmentationType.toString(fragType)+"\n");
		sb.append(" -ptol "+precursorTolerance.getPpmTolerance()+"\n");
		sb.append(" -ftol "+fragmentTolerance.getPpmTolerance()+"\n");
		sb.append(" -enzyme "+enzyme.getName()+"\n");
		sb.append(" -minLength "+minPeptideLength+"\n");
		sb.append(" -maxLength "+maxPeptideLength+"\n");
		sb.append(" -maxMissedCleavage "+maxMissedCleavages+"\n");
		sb.append(" -minCharge "+minCharge+"\n");
		sb.append(" -maxCharge "+maxCharge+"\n");
		sb.append(" -minEluteTime "+minEluteTime+"\n");
		sb.append(" -numberOfReportedPeaks "+numberOfReportedPeaks+"\n");
		sb.append(" -addDecoysToBackground "+addDecoysToBackgound+"\n");
		sb.append(" -dontRunDecoys "+dontRunDecoys+"\n");
		sb.append(" -percolatorThreshold "+percolatorThreshold+"\n");
		
		return sb.toString();
	}
	
	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme, int minPeptideLength,
			int maxPeptideLength, int maxMissedCleavages, byte minCharge, byte maxCharge, int minEluteTime, int numberOfReportedPeaks, boolean addDecoysToBackgound, boolean dontRunDecoys, float percolatorThreshold) {
		this.aaConstants=aaConstants;
		this.fragType=fragType;
		this.precursorTolerance=precursorTolerance;
		this.fragmentTolerance=fragmentTolerance;
		this.enzyme=enzyme;
		this.minPeptideLength=minPeptideLength;
		this.maxPeptideLength=maxPeptideLength;
		this.maxMissedCleavages=maxMissedCleavages;
		this.minCharge=minCharge;
		this.maxCharge=maxCharge;
		this.minEluteTime=minEluteTime;
		this.numberOfReportedPeaks=numberOfReportedPeaks;
		this.addDecoysToBackgound=addDecoysToBackgound;
		this.dontRunDecoys=dontRunDecoys;
		this.percolatorThreshold=percolatorThreshold;
	}
	
	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme,
			int maxMissedCleavages, byte minCharge, byte maxCharge) {
		this.aaConstants=aaConstants;
		this.fragType=fragType;
		this.precursorTolerance=precursorTolerance;
		this.fragmentTolerance=fragmentTolerance;
		this.enzyme=enzyme;
		minPeptideLength=5;
		maxPeptideLength=99;
		this.maxMissedCleavages=maxMissedCleavages;
		this.minCharge=minCharge;
		this.maxCharge=maxCharge;
		minEluteTime=12;
		numberOfReportedPeaks=3;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		percolatorThreshold=0.01f;
	}

	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance fragmentTolerance, MassTolerance precursorTolerance, DigestionEnzyme enzyme) {
		this.aaConstants=aaConstants;
		this.fragType=fragType;
		this.fragmentTolerance=fragmentTolerance;
		this.precursorTolerance=precursorTolerance;
		this.enzyme=enzyme;

		minPeptideLength=5;
		maxPeptideLength=99;
		maxMissedCleavages=1;
		minCharge=2;
		maxCharge=3;
		minEluteTime=12;
		numberOfReportedPeaks=3;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		percolatorThreshold=0.01f;
	}

	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance fragmentTolerance, MassTolerance precursorTolerance, DigestionEnzyme enzyme,
			int maxMissedCleavages) {
		this.aaConstants=aaConstants;
		this.fragType=fragType;
		this.fragmentTolerance=fragmentTolerance;
		this.precursorTolerance=precursorTolerance;
		this.enzyme=enzyme;
		this.maxMissedCleavages=maxMissedCleavages;

		minPeptideLength=5;
		maxPeptideLength=99;
		minCharge=2;
		maxCharge=3;
		minEluteTime=12;
		numberOfReportedPeaks=3;
		addDecoysToBackgound=false;
		dontRunDecoys=false;
		percolatorThreshold=0.01f;
	}

	public AminoAcidConstants getAAConstants() {
		return aaConstants;
	}

	public FragmentationType getFragType() {
		return fragType;
	}

	public MassTolerance getFragmentTolerance() {
		return fragmentTolerance;
	}

	public MassTolerance getPrecursorTolerance() {
		return precursorTolerance;
	}

	public DigestionEnzyme getEnzyme() {
		return enzyme;
	}

	public int getMaxMissedCleavages() {
		return maxMissedCleavages;
	}

	public int getMaxPeptideLength() {
		return maxPeptideLength;
	}

	public int getMinPeptideLength() {
		return minPeptideLength;
	}

	public byte getMaxCharge() {
		return maxCharge;
	}

	public byte getMinCharge() {
		return minCharge;
	}

	public int getMinEluteTime() {
		return minEluteTime;
	}

	public int getNumberOfReportedPeaks() {
		return numberOfReportedPeaks;
	}

	public boolean isAddDecoysToBackgound() {
		return addDecoysToBackgound;
	}
	
	public boolean isDontRunDecoys() {
		return dontRunDecoys;
	}
	
	public float getPercolatorThreshold() {
		return percolatorThreshold;
	}
}
