package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.HashSet;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class AnnotatedLibraryEntry extends LibraryEntry {
	private final FragmentIon[] ionAnnotations;
	private final boolean isDecoy;

	public AnnotatedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime,
			float score, double[] massArray, float[] intensityArray, float[] correlationArray, FragmentIon[] ionAnnotations, boolean isDecoy) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray);
		this.ionAnnotations=ionAnnotations;
		this.isDecoy=isDecoy;
	}

	public AnnotatedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime,
			float score, double[] massArray, float[] intensityArray, float[] correlationArray, FragmentIon[] ionAnnotations) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray);
		this.ionAnnotations=ionAnnotations;
		this.isDecoy=false;
	}

	public AnnotatedLibraryEntry(LibraryEntry entry, SearchParameters parameters) {
		super(entry.getSource(), entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getCopies(),
				entry.getRetentionTime(), entry.getScore(), entry.getMassArray(), entry.getIntensityArray(), entry.getCorrelationArray());

		double[] massArray=entry.getMassArray();
		this.ionAnnotations=new FragmentIon[massArray.length];
		this.isDecoy=false;

		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true)) {
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
		this.isDecoy=false;

		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true)) {
			int[] indicies=parameters.getFragmentTolerance().getIndicies(massArray, fragmentIon.mass);
			for (int i=0; i<indicies.length; i++) {
				ionAnnotations[indicies[i]]=fragmentIon;
			}
		}
	}
	
	public static AnnotatedLibraryEntry getAnnotationsOnly(LibraryEntry entry, SearchParameters parameters) {
		double[] massArray=entry.getMassArray();
		float[] intensityArray=entry.getIntensityArray();
		float[] correlationArray=entry.getCorrelationArray();
		
		TDoubleArrayList newMasses=new TDoubleArrayList();
		TFloatArrayList newIntensities=new TFloatArrayList();
		TFloatArrayList newCorrelations=new TFloatArrayList();
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), true)) {
			int[] indicies=parameters.getFragmentTolerance().getIndicies(massArray, fragmentIon.mass);
			for (int i=0; i<indicies.length; i++) {
				newMasses.add(massArray[indicies[i]]);
				newIntensities.add(intensityArray[indicies[i]]);
				newCorrelations.add(correlationArray[indicies[i]]);
			}
		}
		
		LibraryEntry newEntry=new LibraryEntry(entry.getSource(), entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getCopies(),
				entry.getRetentionTime(), entry.getScore(), newMasses.toArray(), newIntensities.toArray(), newCorrelations.toArray());

		return new AnnotatedLibraryEntry(newEntry, parameters);
	}
	
	@Override
	public boolean isDecoy() {
		return isDecoy;
	}
	/**
	 * 
	 * @return null entries are expected for unannotated peaks!
	 */
	public FragmentIon[] getIonAnnotations() {
		return ionAnnotations;
	}
}
