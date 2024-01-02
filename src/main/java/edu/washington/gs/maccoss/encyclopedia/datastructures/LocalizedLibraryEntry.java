package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.AmbiguousPeptideModSeq;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;

public class LocalizedLibraryEntry extends ChromatogramLibraryEntry {
	private final AmbiguousPeptideModSeq peptideAnnotation;
	private final float localizationScore;
	private final FragmentIon[] localizationIons;
	private final int numberOfModifiableResidues;
	private final int numberOfModifications;
	private final boolean isFullyLocalized;

	public LocalizedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ,
			byte precursorCharge, String peptideModSeq, int copies, float retentionTime, float score,
			double[] massArray, float[] intensityArray, float[] correlationArray, boolean[] quantitativeIonsArray, float[] medianChromatogram,
			Range range, AmbiguousPeptideModSeq peptideAnnotation, float localizationScore, FragmentIon[] localizationIons,
			int numberOfModifiableResidues, int numberOfModifications, boolean isFullyLocalized, Optional<Float> ionMobility, AminoAcidConstants aaConstants) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime,
				score, massArray, intensityArray, correlationArray, quantitativeIonsArray, medianChromatogram, range, ionMobility, aaConstants);
		this.peptideAnnotation = peptideAnnotation;
		this.localizationScore = localizationScore;
		this.localizationIons = localizationIons;
		this.numberOfModifiableResidues = numberOfModifiableResidues;
		this.numberOfModifications = numberOfModifications;
		this.isFullyLocalized = isFullyLocalized;
	}

	public AmbiguousPeptideModSeq getPeptideAnnotation() {
		return peptideAnnotation;
	}

	public float getLocalizationScore() {
		return localizationScore;
	}

	public FragmentIon[] getLocalizationIons() {
		return localizationIons;
	}

	public int getNumberOfModifiableResidues() {
		return numberOfModifiableResidues;
	}

	public int getNumberOfModifications() {
		return numberOfModifications;
	}

	public boolean isFullyLocalized() {
		return isFullyLocalized;
	}
}
