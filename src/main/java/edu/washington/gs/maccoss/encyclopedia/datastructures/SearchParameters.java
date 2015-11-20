package edu.washington.gs.maccoss.encyclopedia.datastructures;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

//@Immutable
public class SearchParameters {
	private final FragmentationType fragType;
	private final MassTolerance precursorTolerance;
	private final MassTolerance fragmentTolerance;
	private final DigestionEnzyme enzyme;
	private final int minPeptideLength=5;
	private final int maxPeptideLength=40;
	private final int maxMissedCleavages=1;
	private final byte minCharge=2;
	private final byte maxCharge=3;
	private final int minEluteTime=12;

	public SearchParameters(FragmentationType fragType, MassTolerance fragmentTolerance, MassTolerance precursorTolerance, DigestionEnzyme enzyme) {
		this.fragType=fragType;
		this.fragmentTolerance=fragmentTolerance;
		this.precursorTolerance=precursorTolerance;
		this.enzyme=enzyme;
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
}
