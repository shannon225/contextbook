package edu.washington.gs.maccoss.encyclopedia.algorithms.phospho;

import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryBackgroundInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;

/**
 * use only for testing
 * @author searleb
 *
 */
public class UnitBackgroundFrequencyCalculator implements BackgroundFrequencyInterface, LibraryBackgroundInterface {
	private final float frequency;

	public UnitBackgroundFrequencyCalculator(float frequency) {
		this.frequency=frequency;
	}
	
	@Override
	public float getFraction(double mass) {
		return frequency;
	}

	@Override
	public Pair<double[], float[]> getRoundedMassCounters(double precursorMz, MassTolerance tolerance) {
		return new Pair<double[], float[]>(new double[] {precursorMz}, new float[] {frequency});
	}

	@Override
	public float[] getFrequencies(double[] ions, double precursorMz, MassTolerance tolerance) {
		float[] frequencies=new float[ions.length];
		Arrays.fill(frequencies, frequency);
		return frequencies;
	}

	@Override
	public LibraryBackgroundInterface getLibraryBackground(double precursorMz, MassTolerance tolerance) {
		return new LibraryBackgroundInterface() {
			@Override
			public float getFraction(double mass) {
				return getFrequencies(new double[] {mass}, precursorMz, tolerance)[0];
			}
		};
	}
	

}
