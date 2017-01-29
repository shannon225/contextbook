package edu.washington.gs.maccoss.encyclopedia.algorithms.xcordia;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Stripe;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.SparseXCorrSpectrum;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.Spectrum;

public class XCorrStripe implements Spectrum {
	private final SparseXCorrSpectrum xcorrSpectrum;
	private final Stripe wrapper;

	public XCorrStripe(SparseXCorrSpectrum xcorrSpectrum, Stripe wrapper) {
		this.xcorrSpectrum=xcorrSpectrum;
		this.wrapper=wrapper;
	}

	@Override
	public String getSpectrumName() {
		return wrapper.getSpectrumName();
	}

	@Override
	public float getScanStartTime() {
		return wrapper.getScanStartTime();
	}

	@Override
	public double getPrecursorMZ() {
		return wrapper.getPrecursorMZ();
	}

	@Override
	public double[] getMassArray() {
		return wrapper.getMassArray();
	}

	@Override
	public float[] getIntensityArray() {
		return wrapper.getIntensityArray();
	}

	@Override
	public float getTIC() {
		return wrapper.getTIC();
	}

	public SparseXCorrSpectrum getXcorrSpectrum() {
		return xcorrSpectrum;
	}
}
