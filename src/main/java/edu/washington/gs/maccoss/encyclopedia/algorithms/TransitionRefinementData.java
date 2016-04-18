package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;

public class TransitionRefinementData {
	private final double[] fragmentMassArray;
	private final ArrayList<float[]> chromatograms;
	private final float[] correlationArray;
	private final float[] integrationArray;
	private final float[] medianChromatogram;
	private final Range range;
	private final Optional<double[]> massArray;
	private final Optional<float[]> intensityArray;
	private final Optional<float[]> rtArray;
	
	
	public TransitionRefinementData(double[] fragmentMassArray, ArrayList<float[]> chromatograms, float[] correlationArray, float[] integrationArray, float[] medianChromatogram, Range range) {
		this(fragmentMassArray, chromatograms, correlationArray, integrationArray, medianChromatogram, range, null, null, null);
	}

	/**
	 * 
	 * @param correlationArray
	 * @param integrationArray
	 * @param medianChromatogram
	 * @param range
	 * @param massArray CAN BE NULL
	 * @param intensityArray CAN BE NULL
	 */
	public TransitionRefinementData(double[] fragmentMassArray, ArrayList<float[]> chromatograms, float[] correlationArray, float[] integrationArray, float[] medianChromatogram, Range range, double[] massArray, float[] intensityArray, float[] rtArray) {
		this.fragmentMassArray=fragmentMassArray;
		this.chromatograms=chromatograms;
		this.correlationArray=correlationArray;
		this.integrationArray=integrationArray;
		this.medianChromatogram=medianChromatogram;
		this.range=range;
		this.massArray=Optional.ofNullable(massArray);
		this.intensityArray=Optional.ofNullable(intensityArray);
		this.rtArray=Optional.ofNullable(rtArray);
	}
	
	public float getTotalIntensity(float minimumCorrelation) {
		float total=0.0f;
		for (int i=0; i<correlationArray.length; i++) {
			if (correlationArray[i]>=minimumCorrelation) {
				total+=integrationArray[i];
			}
		}
		return total;
	}
	public int getTotalQuantIons(float minimumCorrelation) {
		int total=0;
		for (int i=0; i<correlationArray.length; i++) {
			if (correlationArray[i]>=minimumCorrelation) {
				total++;
			}
		}
		return total;
	}
	
	public TransitionRefinementData addPeakData(double[] mass, float[] intensity, float[] rts) {
		return new TransitionRefinementData(fragmentMassArray, chromatograms, correlationArray, integrationArray, medianChromatogram, range, mass, intensity, rts);
	}
	
	public double[] getFragmentMassArray() {
		return fragmentMassArray;
	}
	
	public ArrayList<float[]> getChromatograms() {
		return chromatograms;
	}
	
	public float[] getCorrelationArray() {
		return correlationArray;
	}
	public float[] getIntegrationArray() {
		return integrationArray;
	}
	public float[] getMedianChromatogram() {
		return medianChromatogram;
	}
	public Range getRange() {
		return range;
	}
	public Optional<float[]> getIntensityArray() {
		return intensityArray;
	}
	public Optional<double[]> getMassArray() {
		return massArray;
	}
	public Optional<float[]> getRtArray() {
		return rtArray;
	}
}
