package edu.washington.gs.maccoss.encyclopedia.datastructures;

import java.awt.Color;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.GraphType;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTraceInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class IMSSpectrumWrapper implements Spectrum, XYTraceInterface {
	private final double[] masses;
	private final float[] intensities;
	private final Optional<float[]> ionMobilityArray;
	private final float tic;
	private final float scanStartTime;
	private final String name;
	private final double mz;

	public IMSSpectrumWrapper(Spectrum s) {
		masses=s.getMassArray();
		intensities=s.getIntensityArray();
		ionMobilityArray=s.getIonMobilityArray();
		tic=General.sum(intensities);
		scanStartTime=s.getScanStartTime();
		name=s.getSpectrumName();
		mz=s.getPrecursorMZ();
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
		if (getIonMobilityArray().isPresent()) return GraphType.imsspectrum;
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
