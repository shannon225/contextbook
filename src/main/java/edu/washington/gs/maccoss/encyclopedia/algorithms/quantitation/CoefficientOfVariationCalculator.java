package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.util.ArrayList;
import java.util.HashMap;

import edu.washington.gs.maccoss.encyclopedia.utils.math.General;

public class CoefficientOfVariationCalculator {
	private final HashMap<String, SampleCoordinate> sampleKey;
	private final String[] sampleNames;
	private final float maximumAcceptedCV;

	public CoefficientOfVariationCalculator(HashMap<String, SampleCoordinate> sampleKey, String[] sampleNames, float maximumAcceptedCV) {
		this.sampleKey=sampleKey;
		this.sampleNames=sampleNames;
		this.maximumAcceptedCV=maximumAcceptedCV;
	}
	
	public float getMaximumAcceptedCV() {
		return maximumAcceptedCV;
	}

	public boolean passesCV(ArrayList<String> sourceFiles, float[] intensities) {
		return getCV(sourceFiles, intensities)<=maximumAcceptedCV;
	}

	public float getCV(ArrayList<String> sourceFiles, float[] intensities) {
		float[] sampleAverages=new float[sampleNames.length];
		int[] sampleN=new int[sampleNames.length];
		
		for (int i=0; i<intensities.length; i++) {
			SampleCoordinate coord=sampleKey.get(sourceFiles.get(i));
			if (coord==null) {
				System.err.println("Can't find ["+sourceFiles.get(i)+"]!");
				System.err.println("Keys: {");
				for (String name : sampleKey.keySet()) {
					System.err.println("    "+name);
				}
				System.err.println("}");
			}
			sampleAverages[coord.getSampleIndex()]+=intensities[i];
			sampleN[coord.getSampleIndex()]++;
		}
		
		float[] lmNormalized=new float[intensities.length];
		for (int i=0; i<intensities.length; i++) {
			SampleCoordinate coord=sampleKey.get(sourceFiles.get(i));
			lmNormalized[i]=intensities[i]/sampleAverages[coord.getSampleIndex()];
		}
		
		return General.stdev(lmNormalized);
	}
}
