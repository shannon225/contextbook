package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;

import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class CoefficientOfVariationCalculator {
	private final HashMap<String, SampleCoordinate> sampleKey;
	private final String[] sampleNames;
	private final float maximumAcceptedCV;
	private final TIntObjectHashMap<ArrayList<String>> sampleNamesBySampleIndex;

	public CoefficientOfVariationCalculator(HashMap<String, SampleCoordinate> sampleKey, String[] sampleNames, float maximumAcceptedCV) {
		this.sampleKey=sampleKey;
		this.sampleNames=sampleNames;
		this.maximumAcceptedCV=maximumAcceptedCV;
		
		sampleNamesBySampleIndex=new TIntObjectHashMap<>();
		for (Entry<String, SampleCoordinate> entry : this.sampleKey.entrySet()) {
			int samp=entry.getValue().getSampleIndex();
			ArrayList<String> list=sampleNamesBySampleIndex.get(samp);
			if (list==null) {
				list=new ArrayList<>();
				sampleNamesBySampleIndex.put(samp, list);
			}
			list.add(entry.getKey());
		}
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
	
//	private HashSet<String> alreadySeen=new HashSet<>();
	public float getReplicateNormalizationFactor(String sourceFile, TObjectFloatHashMap<String> ticBySourceFileMap) {
		int samp=sampleKey.get(sourceFile).getSampleIndex();
		
		float target=0.0f;
		float totalTIC=0.0f;
		final ArrayList<String> samplesInRep=sampleNamesBySampleIndex.get(samp);
		for (String name : samplesInRep) {
			float tic=ticBySourceFileMap.get(name);
			if (name.equals(sourceFile)) {
				target=tic;
			}
			totalTIC+=tic;
		}
		float factor=(totalTIC/samplesInRep.size())/target;
//		if (!alreadySeen.contains(sourceFile)) {
//			System.out.println(sourceFile+"\t"+factor);
//			for (String name : samplesInRep) {
//				System.out.println("\t"+name+" --> "+ticBySourceFileMap.get(name));
//			}
//			alreadySeen.add(sourceFile);
//		}
		return factor;
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
