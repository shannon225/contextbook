package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.phospho.PhosphoLocalizationData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.list.array.TFloatArrayList;

public class TransitionRefinementData {
	private final double[] fragmentMassArray; // every considered ion
	private final ArrayList<float[]> chromatograms; // every considered ion
	private final float[] correlationArray; // every considered ion
	private final float[] integrationArray; // every considered ion
	private final float[] backgroundArray; // every considered ion
	private final Optional<float[]> deltaMassArray; // every considered ion
	
	private final float[] medianChromatogram;
	private final Range range;
	
	private final Optional<double[]> massArray; // for every quantified ion
	private final Optional<float[]> intensityArray; // for every quantified ion
	
	private final Optional<float[]> rtArray;
	private Optional<PhosphoLocalizationData> phosphoLocalizationData;
	
	public TransitionRefinementData(double[] fragmentMassArray, ArrayList<float[]> chromatograms, float[] correlationArray, float[] integrationArray, float[] backgroundArray, float[] medianChromatogram, Range range) {
		this(fragmentMassArray, chromatograms, correlationArray, integrationArray, backgroundArray, medianChromatogram, range, null, null, null, null, null);
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
	TransitionRefinementData(double[] fragmentMassArray, ArrayList<float[]> chromatograms, float[] correlationArray, float[] integrationArray, float[] backgroundArray, float[] medianChromatogram, Range range, float[] deltaMassArray, double[] massArray, float[] intensityArray, float[] rtArray, PhosphoLocalizationData phosphoData) {
		this.fragmentMassArray=fragmentMassArray;
		this.chromatograms=chromatograms;
		this.correlationArray=correlationArray;
		this.integrationArray=integrationArray;
		this.backgroundArray=backgroundArray;
		this.medianChromatogram=medianChromatogram;
		this.range=range;
		this.deltaMassArray=Optional.ofNullable(deltaMassArray);
		this.massArray=Optional.ofNullable(massArray);
		this.intensityArray=Optional.ofNullable(intensityArray);
		this.rtArray=Optional.ofNullable(rtArray);
		phosphoLocalizationData=Optional.ofNullable(phosphoData);
	}
	
	public void setPhosphoLocalizationData(Optional<PhosphoLocalizationData> phosphoLocalizationData) {
		this.phosphoLocalizationData=phosphoLocalizationData;
	}
	public Optional<PhosphoLocalizationData> getPhosphoLocalizationData() {
		return phosphoLocalizationData;
	}
	public Pair<Float, Integer> getTopNIntensity(float minimumCorrelation, int n) {
		TFloatArrayList intensities=new TFloatArrayList();
		for (int i=0; i<correlationArray.length; i++) {
			if (correlationArray[i]>=minimumCorrelation) {
				intensities.add(integrationArray[i]);
			}
		}
		intensities.sort();
		
		float total=0.0f;
		int count=1;
		for (int i=intensities.size()-1; i>=0; i--) {
			total+=intensities.get(i);
			if (count>=n) break;
			count++;
		}
		return new Pair<Float, Integer>(total, count);
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
	
	/**
	 * 
	 * @param mass for every quantified ion
	 * @param deltaMass for every considered ion! Must line up with correlationArray
	 * @param intensity for every quantified ion
	 * @param rts used for plotting
	 * @return
	 */
	public TransitionRefinementData addPeakData(float[] deltaMass, double[] mass, float[] intensity, float[] rts) {
		return new TransitionRefinementData(fragmentMassArray, chromatograms, correlationArray, integrationArray, backgroundArray, medianChromatogram, range, deltaMass, mass, intensity, rts, phosphoLocalizationData.isPresent()?phosphoLocalizationData.get():null);
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
	public float[] getBackgroundArray() {
		return backgroundArray;
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
	public Optional<float[]> getDeltaMassArray() {
		return deltaMassArray;
	}
	public Optional<float[]> getRtArray() {
		return rtArray;
	}
}
