package edu.washington.gs.maccoss.encyclopedia.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.PeptidePrecursor;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.list.array.TFloatArrayList;

public class TransitionRefinementData implements PeptidePrecursor {
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
	
	private final String peptideModSeq;
	private final byte precursorCharge;
	
	// NOTE: with mods, the cannonical form points to the modificationQuantData
	// map of all possible forms. Individual forms are specified as
	// localizationData. Either set one or the other of these! Don't set both!
	private Optional<ModificationLocalizationData> localizationData;
	private Optional<HashMap<String, TransitionRefinementData>> modificationQuantData;
	
	public TransitionRefinementData(String peptideModSeq, byte precursorCharge, double[] fragmentMassArray, ArrayList<float[]> chromatograms, float[] correlationArray, float[] integrationArray, float[] backgroundArray, float[] medianChromatogram, Range range) {
		this(peptideModSeq, precursorCharge, fragmentMassArray, chromatograms, correlationArray, integrationArray, backgroundArray, medianChromatogram, range, null, null, null, null, null, null);
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
	TransitionRefinementData(String peptideModSeq, byte precursorCharge, double[] fragmentMassArray, ArrayList<float[]> chromatograms, float[] correlationArray, float[] integrationArray, float[] backgroundArray, float[] medianChromatogram, Range range, float[] deltaMassArray, double[] massArray, float[] intensityArray, float[] rtArray, ModificationLocalizationData localizationData, HashMap<String, TransitionRefinementData> modificationQuantData) {
		this.peptideModSeq=peptideModSeq;
		this.precursorCharge=precursorCharge;
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
		this.localizationData=Optional.ofNullable(localizationData);
		this.modificationQuantData=Optional.ofNullable(modificationQuantData);
	}
	
	@Override
	public int compareTo(PeptidePrecursor o) {
		if (o==null) return 1;
		int c=getPeptideModSeq().compareTo(o.getPeptideModSeq());
		if (c!=0) return c;
		return Byte.compare(getPrecursorCharge(), o.getPrecursorCharge());
	}
	
	public String getPeptideModSeq() {
		return peptideModSeq;
	}
	
	public String getPeptideSeq() {
		StringBuilder sb=new StringBuilder();
		for (char c : peptideModSeq.toCharArray()) {
			if (Character.isLetter(c)) {
				sb.append(c);
			}
		}
		return sb.toString();
	}
	
	public byte getPrecursorCharge() {
		return precursorCharge;
	}
	
	public void setModificationLocalizationData(Optional<ModificationLocalizationData> newLocalizationData) {
		this.localizationData=newLocalizationData;
	}
	public void setModificationQuantData(Optional<HashMap<String, TransitionRefinementData>> newQuantData) {
		this.modificationQuantData=newQuantData;
	}
	public Optional<ModificationLocalizationData> getLocalizationData() {
		return localizationData;
	}
	public Optional<HashMap<String, TransitionRefinementData>> getModificationQuantData() {
		return modificationQuantData;
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
	
	public float getApexRT() {
		if (!rtArray.isPresent()) {
			throw new EncyclopediaException("Requesting apex RT but no retention times are loaded!");
		}
		
		int bestIndex=-1;
		float bestMedian=-Float.MAX_VALUE;
		for (int i=0; i<medianChromatogram.length; i++) {
			if (medianChromatogram[i]>bestMedian) {
				bestIndex=i;
				bestMedian=medianChromatogram[i];
			}
		}
		return rtArray.get()[bestIndex];
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
		return new TransitionRefinementData(peptideModSeq, precursorCharge, fragmentMassArray, chromatograms, correlationArray, integrationArray, backgroundArray, medianChromatogram, range, deltaMass, mass, intensity, rts, localizationData.isPresent()?localizationData.get():null, modificationQuantData.isPresent()?modificationQuantData.get():null);
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
