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
	
	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme, int minPeptideLength,
			int maxPeptideLength, int maxMissedCleavages, byte minCharge, byte maxCharge, int minEluteTime, int numberOfReportedPeaks, boolean addDecoysToBackgound) {
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
}
