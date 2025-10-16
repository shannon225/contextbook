package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.awt.Color;
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
	private final double[] masses;
	private final float[] intensities;
	private final Optional<float[]> ionMobilityArray;
	private final FragmentIon[] annotations;
	private final float tic;
	private final float scanStartTime;
	private final String name;
	private final double mz;
	
	public AnnotatedSpectrum(String name, double mz, float scanStartTime, double[] masses, float[] intensities, Optional<float[]> ionMobilityArray, FragmentIon[] annotations) {
		super();
		this.masses = masses;
		this.intensities = intensities;
		this.ionMobilityArray=ionMobilityArray;
		this.tic = General.sum(intensities);
		this.scanStartTime = scanStartTime;
		this.name = name;
		this.mz = mz;
		this.annotations=annotations;
	}

	public AnnotatedSpectrum(Spectrum s, PeptidePrecursor entry, SearchParameters parameters) {
		masses=s.getMassArray();
		intensities=s.getIntensityArray();
		ionMobilityArray=s.getIonMobilityArray();
		tic=General.sum(intensities);
		scanStartTime=s.getScanStartTime();
		name=s.getSpectrumName();
		mz=s.getPrecursorMZ();
		
		double[] massArray=s.getMassArray();
		this.annotations=new FragmentIon[massArray.length];

		FragmentationModel model=PeptideUtils.getPeptideModel(entry.getPeptideModSeq(), parameters.getAAConstants());
		for (FragmentIon fragmentIon : model.getPrimaryIonObjects(parameters.getFragType(), entry.getPrecursorCharge(), false)) {
			int[] indicies=parameters.getFragmentTolerance().getIndicies(massArray, fragmentIon.getMass());
			for (int i=0; i<indicies.length; i++) {
				annotations[indicies[i]]=fragmentIon;
			}
		}
	}
	
	public AnnotatedSpectrum(Spectrum s, TDoubleObjectHashMap<String> annotationMap, MassTolerance tolerance) {
		masses=s.getMassArray();
		intensities=s.getIntensityArray();
		ionMobilityArray=s.getIonMobilityArray();
		tic=General.sum(intensities);
		scanStartTime=s.getScanStartTime();
		name=s.getSpectrumName();
		mz=s.getPrecursorMZ();
		
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
