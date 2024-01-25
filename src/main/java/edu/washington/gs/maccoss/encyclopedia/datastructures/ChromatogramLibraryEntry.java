package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class ChromatogramLibraryEntry extends LibraryEntry implements Chromatogram {
	private final float[] medianChromatogram;
	private final Range range;
	
	public ChromatogramLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score, double[] massArray,
			float[] intensityArray, float[] correlationArray, boolean[] quantifiedIonsArray, float[] medianChromatogram, Range range, Optional<Float> ionMobility, AminoAcidConstants aaConstants) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray, quantifiedIonsArray, ionMobility, aaConstants);
		this.medianChromatogram=medianChromatogram;
		this.range=range;
	}
	
	public ChromatogramLibraryEntry(String source, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, 
			String massCorrectedPeptideModSeq, int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, float[] correlationArray, 
			boolean[] quantifiedIonsArray, boolean keepNegativeIntensities, float[] medianChromatogram, Range range, Optional<Float> ionMobility) {
		super(source, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, massCorrectedPeptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray, quantifiedIonsArray, ionMobility, false);
		this.medianChromatogram=medianChromatogram;
		this.range=range;
	}
	
	/**
	 * only use for testing
	 * @param newMassArray
	 * @param newIntensityArray
	 * @param newCorrelationArray
	 * @return
	 */
	@Override
	public ChromatogramLibraryEntry updateMS2(double[] newMassArray, float[] newIntensityArray, float[] newCorrelationArray, boolean[] newQuantifiedIonsArray) {
		return new ChromatogramLibraryEntry(getSource(), getAccessions(), getSpectrumIndex(), getPrecursorMZ(), getPrecursorCharge(), getPeptideModSeq(), 
				getPeptideModSeq(), getCopies(), getRetentionTime(), getScore(), newMassArray, newIntensityArray, newCorrelationArray, 
				newQuantifiedIonsArray, false, getMedianChromatogram(), getRtRange(), getIonMobility());
	}

	public LibraryEntry sqrt() {
		float[] sqrt=General.protectedSqrt(getIntensityArray());
		return new ChromatogramLibraryEntry(getSource(), getAccessions(), getSpectrumIndex(), getPrecursorMZ(), getPrecursorCharge(), getPeptideModSeq(), 
				getPeptideModSeq(), getCopies(), getRetentionTime(), getScore(), getMassArray(), sqrt, getCorrelationArray(), 
				getQuantifiedIonsArray(), false, getMedianChromatogram(), getRtRange(), getIonMobility());
	}
	
	public Range getRtRange() {
		return range;
	}
	
	public float[] getMedianChromatogram() {
		return medianChromatogram;
	}
	
	public float getDurationInSec() {
		return range.getRange();
	}
}
