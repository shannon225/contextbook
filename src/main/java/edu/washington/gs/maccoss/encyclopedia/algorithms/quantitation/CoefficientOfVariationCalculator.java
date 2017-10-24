package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;

public class CoefficientOfVariationCalculator {
	private final HashMap<String, SampleCoordinate> sampleKey;
	private final String[] sampleNames;
	private final float maximumAcceptedCV;

	public CoefficientOfVariationCalculator(HashMap<String, SampleCoordinate> sampleKey, String[] sampleNames, float maximumAcceptedCV) {
		this.sampleKey=sampleKey;
		this.sampleNames=sampleNames;
		this.maximumAcceptedCV=maximumAcceptedCV;
	}
	
	public ArrayList<String> getSortedSampleNames() {
		TreeMap<SampleCoordinate, String> samples=new TreeMap<>();
		for (Entry<String, SampleCoordinate> entry : sampleKey.entrySet()) {
			samples.put(entry.getValue(), entry.getKey());
		}
		return new ArrayList<String>(samples.values());
	}
	
	public float getMaximumAcceptedCV() {
		return maximumAcceptedCV;
	}

	/**
	 * 
	 * @param sourceFiles
	 * @param intensities
	 * @return {CV, atLeastSampleFullyMeasured}
	 */
	public Pair<Float, Boolean> getCV(ArrayList<String> sourceFiles, float[] intensities) {
		float[] sampleAverages=new float[sampleNames.length];
		int[] sampleN=new int[sampleNames.length];
		boolean[] isFullyMeasured=new boolean[sampleAverages.length];
		Arrays.fill(isFullyMeasured, true);
		
		for (int i=0; i<intensities.length; i++) {
			SampleCoordinate coord=sampleKey.get(sourceFiles.get(i));
			if (coord==null) {
				Logger.errorLine("Can't find ["+sourceFiles.get(i)+"]!");
				Logger.errorLine("Keys: {");
				for (String name : sampleKey.keySet()) {
					Logger.errorLine("    "+name);
				}
				Logger.errorLine("}");
			}
			sampleAverages[coord.getSampleIndex()]+=intensities[i];
			sampleN[coord.getSampleIndex()]++;
			if (intensities[i]==0.0f) {
				isFullyMeasured[coord.getSampleIndex()]=false;
			}
		}
		
		TFloatArrayList lmNormalized=new TFloatArrayList();
		for (int i=0; i<intensities.length; i++) {
			SampleCoordinate coord=sampleKey.get(sourceFiles.get(i));
			if (sampleAverages[coord.getSampleIndex()]>0) {
				lmNormalized.add(intensities[i]/sampleAverages[coord.getSampleIndex()]);
			}
		}
		
		boolean atLeastSampleFullyMeasured=false;
		for (boolean b : isFullyMeasured) {
			if (b) {
				atLeastSampleFullyMeasured=true;
				break;
			}
		}
		
		return new Pair<Float, Boolean>(General.stdev(lmNormalized.toArray()), atLeastSampleFullyMeasured);
	}
}
