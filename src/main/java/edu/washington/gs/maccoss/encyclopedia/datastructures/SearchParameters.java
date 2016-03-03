package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.DigestionEnzyme;

public class SearchParameters {

	protected final AminoAcidConstants aaConstants;
	protected final FragmentationType fragType;
	protected final MassTolerance precursorTolerance;
	protected final MassTolerance fragmentTolerance;
	protected final DigestionEnzyme enzyme;
	protected final float percolatorThreshold;
	protected final File percolatorLocation;
	protected final boolean deconvoluteOverlappingWindows;
	protected final int numberOfThreadsUsed;	
	protected final float targetWindowCenter;
	protected final float expectedPeakWidth;
	protected final boolean runPhosphoLocalization;

	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme,
			float percolatorThreshold, File percolatorLocation, boolean deconvoluteOverlappingWindows, int numberOfThreadsUsed, float expectedPeakWidth, float targetWindowCenter, boolean runPhosphoLocalization) {
		this.aaConstants=aaConstants;
		this.fragType=fragType;
		this.precursorTolerance=precursorTolerance;
		this.fragmentTolerance=fragmentTolerance;
		this.enzyme=enzyme;
		this.percolatorThreshold=percolatorThreshold;
		this.percolatorLocation=percolatorLocation;
		this.deconvoluteOverlappingWindows=deconvoluteOverlappingWindows;
		this.numberOfThreadsUsed=numberOfThreadsUsed;
		this.expectedPeakWidth=expectedPeakWidth;
		this.targetWindowCenter=targetWindowCenter;
		this.runPhosphoLocalization=runPhosphoLocalization;
	}

	public String toString() {
		final StringBuilder sb=new StringBuilder();
		sb.append(" -fixed "+aaConstants.getFixedModString()+"\n");
		sb.append(" -frag "+FragmentationType.toString(fragType)+"\n");
		sb.append(" -ptol "+precursorTolerance.getPpmTolerance()+"\n");
		sb.append(" -ftol "+fragmentTolerance.getPpmTolerance()+"\n");
		sb.append(" -enzyme "+enzyme.getName()+"\n");
		sb.append(" -percolatorThreshold "+percolatorThreshold+"\n");
		sb.append(" -percolatorLocation "+percolatorLocation+"\n");
		sb.append(" -deconvoluteOverlappingWindows "+deconvoluteOverlappingWindows+"\n");
		sb.append(" -numberOfThreadsUsed "+numberOfThreadsUsed+"\n");
		sb.append(" -expectedPeakWidth "+expectedPeakWidth+"\n");
		sb.append(" -runPhosphoLocalization "+runPhosphoLocalization+"\n");
		if (useTargetWindowCenter()) {
			sb.append(" -targetWindowCenter "+targetWindowCenter+"\n");
		}
		return sb.toString();
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

	public float getPercolatorThreshold() {
		return percolatorThreshold;
	}

	public Optional<File> getPercolatorLocation() {
		return Optional.ofNullable(percolatorLocation);
	}

	public boolean isDeconvoluteOverlappingWindows() {
		return deconvoluteOverlappingWindows;
	}
	public int getNumberOfThreadsUsed() {
		return numberOfThreadsUsed;
	}
	public float getTargetWindowCenter() {
		return targetWindowCenter;
	}
	public boolean useTargetWindowCenter() {
		return targetWindowCenter>0;
	}
	
	public float getExpectedPeakWidth() {
		return expectedPeakWidth;
	}
	
	public boolean isRunPhosphoLocalization() {
		return runPhosphoLocalization;
	}
}