package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class AnnotatedLibraryEntry extends LibraryEntry {
	private final FragmentIon[] ionAnnotations;

	public AnnotatedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime,
			float score, double[] massArray, float[] intensityArray, float[] correlationArray, FragmentIon[] ionAnnotations) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray);
		this.ionAnnotations=ionAnnotations;
	}

	public AnnotatedLibraryEntry(LibraryEntry entry, SearchParameters parameters) {
		super(entry.getSource(), entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getCopies(),
				entry.getRetentionTime(), entry.getScore(), entry.getMassArray(), entry.getIntensityArray(), entry.getCorrelationArray());

		double[] massArray=entry.getMassArray();
		this.ionAnnotations=new FragmentIon[massArray.length];

		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge())) {
			int[] indicies=parameters.getFragmentTolerance().getIndicies(massArray, fragmentIon.mass);
			for (int i=0; i<indicies.length; i++) {
				ionAnnotations[indicies[i]]=fragmentIon;
			}
		}
	}

	public AnnotatedLibraryEntry(PeptidePrecursor entry, Spectrum spectrum, SearchParameters parameters) {
		super(spectrum.getSpectrumName(), new HashSet<String>(), 1, parameters.getAAConstants().getChargedMass(entry.getPeptideModSeq(), entry.getPrecursorCharge()), entry.getPrecursorCharge(),
				entry.getPeptideModSeq(), 1, spectrum.getScanStartTime(), 0.0f, spectrum.getMassArray(), spectrum.getIntensityArray(), new float[spectrum.getMassArray().length]);

		double[] massArray=spectrum.getMassArray();
		this.ionAnnotations=new FragmentIon[massArray.length];

		FragmentationModel model=new FragmentationModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge())) {
			int[] indicies=parameters.getFragmentTolerance().getIndicies(massArray, fragmentIon.mass);
			for (int i=0; i<indicies.length; i++) {
				ionAnnotations[indicies[i]]=fragmentIon;
			}
		}
	}
	

	/**
	 * 
	 * @return null entries are expected for unannotated peaks!
	 */
	public FragmentIon[] getIonAnnotations() {
		return ionAnnotations;
	}
}
