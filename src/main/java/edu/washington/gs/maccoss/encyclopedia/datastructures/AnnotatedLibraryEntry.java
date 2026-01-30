package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Ion;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.ScoredObject;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class AnnotatedLibraryEntry extends LibraryEntry {
	private final FragmentIon[] ionAnnotations;
	private final boolean isDecoy;

	public AnnotatedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime,
			float score, double[] massArray, float[] intensityArray, float[] correlationArray, boolean[] quantifiedIonsArray, FragmentIon[] ionAnnotations, Optional<Float> ionMobility, boolean isDecoy, AminoAcidConstants aaConstants) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray, quantifiedIonsArray, ionMobility, aaConstants);
		this.ionAnnotations=ionAnnotations;
		this.isDecoy=isDecoy;
	}

	public AnnotatedLibraryEntry(String sourceFile, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, int copies, float retentionTime,
			float score, double[] massArray, float[] intensityArray, float[] correlationArray, boolean[] quantifiedIonsArray, FragmentIon[] ionAnnotations, Optional<Float> ionMobility, AminoAcidConstants aaConstants) {
		super(sourceFile, accessions, spectrumIndex, precursorMZ, precursorCharge, peptideModSeq, copies, retentionTime, score, massArray, intensityArray, correlationArray, quantifiedIonsArray, ionMobility, aaConstants);
		this.ionAnnotations=ionAnnotations;
		this.isDecoy=false;
	}

	public AnnotatedLibraryEntry(LibraryEntry entry, SearchParameters parameters) {
		this(entry, parameters, false);
	}
	public AnnotatedLibraryEntry(LibraryEntry entry, SearchParameters parameters, boolean keepNegativeIntensities) {
		super(entry.getSource(), entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getCopies(),
				entry.getRetentionTime(), entry.getScore(), entry.getMassArray(), entry.getIntensityArray(), entry.getCorrelationArray(), entry.getQuantifiedIonsArray(), entry.getIonMobility(), parameters.getAAConstants(), keepNegativeIntensities);

		double[] massArray=entry.getMassArray();
		this.ionAnnotations=new FragmentIon[massArray.length];
		this.isDecoy=false;

		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false)) {
			int[] indicies=parameters.getFragmentTolerance().getIndices(massArray, fragmentIon.getMass());
			for (int i=0; i<indicies.length; i++) {
				ionAnnotations[indicies[i]]=fragmentIon;
			}
		}
	}

	public AnnotatedLibraryEntry(PeptidePrecursor entry, Spectrum spectrum, SearchParameters parameters) {
		super(spectrum.getSpectrumName(), new HashSet<String>(), 1, parameters.getAAConstants().getChargedMass(entry.getPeptideModSeq(), entry.getPrecursorCharge()), entry.getPrecursorCharge(),
				entry.getPeptideModSeq(), 1, spectrum.getScanStartTime(), 0.0f, spectrum.getMassArray(), spectrum.getIntensityArray(), new float[spectrum.getMassArray().length], new boolean[spectrum.getMassArray().length], LibraryEntry.getAverageIonMobilityFromArray(spectrum.getIonMobilityArray()), parameters.getAAConstants());

		double[] massArray=spectrum.getMassArray();
		this.ionAnnotations=new FragmentIon[massArray.length];
		this.isDecoy=false;

		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false)) {
			int[] indicies=parameters.getFragmentTolerance().getIndices(massArray, fragmentIon.getMass());
			for (int i=0; i<indicies.length; i++) {
				ionAnnotations[indicies[i]]=fragmentIon;
			}
		}
	}
	
	public static AnnotatedLibraryEntry getAnnotationsOnly(LibraryEntry entry, SearchParameters parameters) {
		return getAnnotationsOnly(entry, parameters, 0);
	}

	public static AnnotatedLibraryEntry getAnnotationsOnly(LibraryEntry entry, SearchParameters parameters, int minIonNumber) {
		double[] massArray=entry.getMassArray();
		float[] intensityArray=entry.getIntensityArray();
		float[] correlationArray=entry.getCorrelationArray();
		boolean[] quantifiedIonsArray=entry.getQuantifiedIonsArray();
		
		TDoubleArrayList newMasses=new TDoubleArrayList();
		TFloatArrayList newIntensities=new TFloatArrayList();
		TFloatArrayList newCorrelations=new TFloatArrayList();
		ArrayList<Boolean> newQuantifiedIons=new ArrayList<Boolean>();
		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (Ion fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false)) {
			if (minIonNumber==0||fragmentIon.getIndex()>=minIonNumber) {
				int[] indicies=parameters.getFragmentTolerance().getIndices(massArray, fragmentIon.getMass());
				for (int i=0; i<indicies.length; i++) {
					newMasses.add(massArray[indicies[i]]);
					newIntensities.add(intensityArray[indicies[i]]);
					newCorrelations.add(correlationArray[indicies[i]]);
					newQuantifiedIons.add(quantifiedIonsArray[indicies[i]]);
				}
			}
		}
		boolean[] newQuantifiedIonsArray=new boolean[newQuantifiedIons.size()];
		for (int i = 0; i < newQuantifiedIonsArray.length; i++) {
			newQuantifiedIonsArray[i]=newQuantifiedIons.get(i);
		}
		
		LibraryEntry newEntry=new LibraryEntry(entry.getSource(), entry.getAccessions(), entry.getSpectrumIndex(), entry.getPrecursorMZ(), entry.getPrecursorCharge(), entry.getPeptideModSeq(), entry.getCopies(),
				entry.getRetentionTime(), entry.getScore(), newMasses.toArray(), newIntensities.toArray(), newCorrelations.toArray(), newQuantifiedIonsArray, entry.getIonMobility(), parameters.getAAConstants());

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
	
	/**
	 * 
	 * @param n
	 * @param minIndex
	 * @param minimumPercent
	 * @param mask can be null
	 * @return
	 */
	public FragmentIon[] getMostCorrelatedAnnotatedIons(int n, int minIndex, float minimumCorrelation) {
		float[] correlationArray=getCorrelationArray();
		ArrayList<ScoredObject<FragmentIon>> scored=new ArrayList<ScoredObject<FragmentIon>>();
		for (int i=0; i<ionAnnotations.length; i++) {
			if (ionAnnotations[i]!=null) {
				if (ionAnnotations[i].getIndex()>=minIndex&&correlationArray[i]>minimumCorrelation) {
					scored.add(new ScoredObject<FragmentIon>(correlationArray[i], ionAnnotations[i]));
				}
			}
		}
		scored.sort(null);
		
		ArrayList<FragmentIon> best=new ArrayList<FragmentIon>();
		for (int i=scored.size()-1; i>=0; i--) {
			best.add(scored.get(i).y);
			if (best.size()>=n) break;
		}
		return best.toArray(new FragmentIon[0]);
	}
	
	/**
	 * 
	 * @param n
	 * @param minIndex
	 * @param minimumPercent
	 * @param mask can be null
	 * @return
	 */
	public FragmentIon[] getMostIntenseAnnotatedIons(int n, int minIndex, float minimumPercent, boolean[] mask) {
		float[] intensityArray=getIntensityArray();
		ArrayList<ScoredObject<FragmentIon>> scored=new ArrayList<ScoredObject<FragmentIon>>();
		float maxIntensity=0.0f;
		for (int i=0; i<ionAnnotations.length; i++) {
			if (ionAnnotations[i]!=null&&(mask==null||mask[i])) {
				if (ionAnnotations[i].getIndex()>=minIndex) {
					if (maxIntensity<intensityArray[i]) {
						maxIntensity=intensityArray[i];
					}
					scored.add(new ScoredObject<FragmentIon>(intensityArray[i], ionAnnotations[i]));
				}
			}
		}
		scored.sort(null);
		
		ArrayList<FragmentIon> best=new ArrayList<FragmentIon>();
		for (int i=scored.size()-1; i>=0; i--) {
			if (scored.get(i).x/maxIntensity<minimumPercent) break;
			best.add(scored.get(i).y);
			if (best.size()>=n) break;
		}
		return best.toArray(new FragmentIon[0]);
	}
}
