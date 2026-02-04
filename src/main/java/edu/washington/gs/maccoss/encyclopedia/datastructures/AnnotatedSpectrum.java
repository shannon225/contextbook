package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.awt.Color;
import java.util.HashSet;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.AnnotatedFragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.FragmentIon;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import gnu.trove.procedure.TDoubleObjectProcedure;

public class AnnotatedSpectrum implements Spectrum, XYTraceInterface {
	private final String peptideModSeq;
	private final double[] masses;
	private final float[] intensities;
	private final Optional<float[]> ionMobilityArray;
	private final FragmentIon[] annotations;
	private final float tic;
	private final float scanStartTime;
	private final String name;
	private final double mz;
	private final byte precursorCharge;
	
	public LibraryEntry getEntry() {
		double precursorMZ=new AminoAcidConstants().getChargedMass(peptideModSeq, precursorCharge);
		//String source, HashSet<String> accessions, int spectrumIndex, double precursorMZ, byte precursorCharge, String peptideModSeq, String massCorrectedPeptideModSeq, 
		//int copies, float retentionTime, float score, double[] massArray, float[] intensityArray, float[] correlationArray, boolean[] quantifiedIonsArray, Optional<Float> ionMobility
		Optional<Float> ionMobility;
		if (ionMobilityArray.isPresent()) {
			ionMobility=Optional.of(Float.valueOf(General.sum(General.multiply(intensities, ionMobilityArray.get()))/tic));
		} else {
			ionMobility=Optional.empty();
		}
		return new LibraryEntry(name, new HashSet<String>(), 1, precursorMZ, precursorCharge, peptideModSeq, peptideModSeq, 
				1, scanStartTime, 0.0f, masses, intensities, new float[masses.length], General.getBooleanUnitArray(masses.length, true), ionMobility, false);
	}

	public AnnotatedSpectrum(Spectrum s, PeptidePrecursor entry, SearchParameters parameters) {
		peptideModSeq=entry.getPeptideModSeq();
		masses=s.getMassArray();
		intensities=s.getIntensityArray();
		ionMobilityArray=s.getIonMobilityArray();
		tic=General.sum(intensities);
		scanStartTime=s.getScanStartTime();
		name=s.getSpectrumName();
		mz=s.getPrecursorMZ();
		precursorCharge=entry.getPrecursorCharge();
		
		double[] massArray=s.getMassArray();
		this.annotations=new FragmentIon[massArray.length];

		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false)) {
			int[] indicies=parameters.getFragmentTolerance().getIndices(massArray, fragmentIon.getMass());
			for (int i=0; i<indicies.length; i++) {
				annotations[indicies[i]]=fragmentIon;
			}
		}
	}
	
	public AnnotatedSpectrum(String peptideModSeq, Spectrum s, TDoubleObjectHashMap<String> annotationMap, byte precursorCharge, MassTolerance tolerance) {
		this.peptideModSeq=peptideModSeq;
		masses=s.getMassArray();
		intensities=s.getIntensityArray();
		ionMobilityArray=s.getIonMobilityArray();
		tic=General.sum(intensities);
		scanStartTime=s.getScanStartTime();
		name=s.getSpectrumName();
		mz=s.getPrecursorMZ();
		this.precursorCharge=precursorCharge;
		
		annotations=new FragmentIon[masses.length];
		annotationMap.forEachEntry(new TDoubleObjectProcedure<String>() {
			@Override
			public boolean execute(double a, String b) {
				Optional<Integer> index=tolerance.getIndex(masses, a);
				if (index.isPresent()) {
					annotations[index.get()]=new AnnotatedFragmentIon(a, b);
				}
				return true;
			}
		});
	}
	
	@Override
	public Optional<Color> getColor() {
		return Optional.ofNullable((Color)null);
	}
	@Override
	public String getName() {
		return name;
	}
	@Override
	public Optional<Float> getThickness() {
		return Optional.ofNullable((Float)null);
	}
	@Override
	public GraphType getType() {
		return GraphType.spectrum;
	}
	@Override
	public Pair<double[], double[]> toArrays() {
		return new Pair<double[], double[]>(masses, General.toDoubleArray(intensities));
	}
	
	@Override
	public int size() {
		return masses.length;
	}
	
	public FragmentIon[] getAnnotations() {
		return annotations;
	}
	
	public Optional<float[]> getIonMobilityArray() {
		return ionMobilityArray;
	}
	
	public double[] getMassArray() {
		return masses;
	}
	public float[] getIntensityArray() {
		return intensities;
	}
	public float getTIC() {
		return tic;
	}
	public float getScanStartTime() {
		return scanStartTime;
	}
	public String getSpectrumName() {
		return name;
	}
	public double getPrecursorMZ() {
		return mz;
	}
}
