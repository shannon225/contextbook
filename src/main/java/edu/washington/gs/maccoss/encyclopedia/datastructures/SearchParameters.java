package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.io.File;

import com.google.common.base.Optional;

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

	public SearchParameters(AminoAcidConstants aaConstants, FragmentationType fragType, MassTolerance precursorTolerance, MassTolerance fragmentTolerance, DigestionEnzyme enzyme,
			float percolatorThreshold, File percolatorLocation, boolean deconvoluteOverlappingWindows) {
		this.aaConstants=aaConstants;
		this.fragType=fragType;
		this.precursorTolerance=precursorTolerance;
		this.fragmentTolerance=fragmentTolerance;
		this.enzyme=enzyme;
		this.percolatorThreshold=percolatorThreshold;
		this.percolatorLocation=percolatorLocation;
		this.deconvoluteOverlappingWindows=deconvoluteOverlappingWindows;
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
		return Optional.fromNullable(percolatorLocation);
	}

	public boolean isDeconvoluteOverlappingWindows() {
		return deconvoluteOverlappingWindows;
	}

}