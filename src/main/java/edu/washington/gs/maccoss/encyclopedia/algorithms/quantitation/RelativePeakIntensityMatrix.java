package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYTrace;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.XYZPoint;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.MassTolerance;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeakFrequencyCalculator;
import edu.washington.gs.maccoss.encyclopedia.utils.math.General;
import edu.washington.gs.maccoss.encyclopedia.utils.math.QuickMedian;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class RelativePeakIntensityMatrix {
	private static final float MAXIMUM_ALLOWED_SUM_SQUARED_ERRORS = 0.05f;
	private final MassTolerance tolerance;
	private final PeakFrequencyCalculator fragmentList;
	private final Optional<LibraryEntry> target;
	private final HashMap<String, LibraryEntry> entryMap=new HashMap<>();
	
	
	public RelativePeakIntensityMatrix(MassTolerance tolerance, Optional<LibraryEntry> target) {
		this.tolerance=tolerance;
		this.fragmentList=new PeakFrequencyCalculator(tolerance);
		this.target=target;
	}

	public void addPeak(String sample, LibraryEntry entry) {
		entryMap.put(sample, entry);
		double[] masses=entry.getMassArray();
		float[] intensity=entry.getIntensityArray();
		float[] correlation=entry.getCorrelationArray();
		for (int i=0; i<correlation.length; i++) {
			if (correlation[i]>=TransitionRefiner.quantitativeCorrelationThreshold) {
				fragmentList.increment(masses[i], intensity[i]);
			}
		}
	}
	
	/**
	 * 
	 * @param keys must be the same as "sample" in addPeak
	 * @return
	 */
	public float[] integratePeptides(double[] nBestPeaks, ArrayList<String> keys) {
		float[] results=new float[keys.size()];
		for (int i = 0; i < results.length; i++) {
			LibraryEntry entry=entryMap.get(keys.get(i));
			if (entry!=null) {
				results[i]=General.sum(tolerance.getIntegratedIntensities(entry.getMassArray(), entry.getIntensityArray(), nBestPeaks));
			}
		}
		return results;
	}
	
	public double[] pickNBestPeaks(int targetNumberOfPoints, int minNumberOfPoints) {
		ArrayList<String> keys=new ArrayList<>(entryMap.keySet());
		Collections.sort(keys);
		
		// check edge cases for 0 and 1 samples
		Optional<double[]> edgeCheck=checkEdgeCases(targetNumberOfPoints);
		if (edgeCheck.isPresent()) return edgeCheck.get();
		
		// get list of fragments that were seen in at least one sample above the quantitativeCorrelationThreshold
		double[] possibleMasses=XYTrace.toArrays(fragmentList.toPoints()).x;
		Arrays.sort(possibleMasses);
		
		// check edge case for no ions
		if (possibleMasses.length==0) return new double[0]; 
		if (keys.size()==0) return new double[0];

		float[][] prelimDataMatrix=generateDataMatrix(possibleMasses, keys); // sample X fragment
		float[] targetIntensities;
		if (target.isPresent()) {
			targetIntensities=getPossibleIntensities(possibleMasses, target.get());
		} else {
			// generate data matrix of all intensities for those fragments
			float[][] prelimTransposedDataMatrix=General.transposeMatrix(prelimDataMatrix); // fragment X sample
			float[] prelimMedianRelativeIntensities=new float[prelimTransposedDataMatrix.length];
			for (int i = 0; i < prelimMedianRelativeIntensities.length; i++) {
				prelimMedianRelativeIntensities[i]=QuickMedian.median(prelimTransposedDataMatrix[i]);
			}
			targetIntensities=prelimMedianRelativeIntensities;
		}

		// remove scored peptides that do not match the median fragment intensity pattern
		float[] sampleSumSquaredErrors=new float[keys.size()];
		ArrayList<String> badSamples=new ArrayList<>();
		for (int i = 0; i < keys.size(); i++) {
			float[] thisSamplesIntensities=prelimDataMatrix[i];
			
			for (int j = 0; j < thisSamplesIntensities.length; j++) {
				float delta=thisSamplesIntensities[j]-targetIntensities[j];
				sampleSumSquaredErrors[i]+=delta*delta;
			}
			if (sampleSumSquaredErrors[i]>MAXIMUM_ALLOWED_SUM_SQUARED_ERRORS) {
				badSamples.add(keys.get(i));
			}
		}
		
		if (keys.size()==badSamples.size()) {
			// find the best sample and use that
			float minimumSSE=Float.MAX_VALUE;
			int index=0;
			for (int i = 0; i < sampleSumSquaredErrors.length; i++) {
				if (sampleSumSquaredErrors[i]<minimumSSE) {
					minimumSSE=sampleSumSquaredErrors[i];
					index=i;
				}
			}
			badSamples.remove(keys.get(index));
		}
		
		for (String key : badSamples) {
			keys.remove(key);
		}
		
		// regenerate matrix with only the "good" peptides
		float[][] dataMatrix=generateDataMatrix(possibleMasses, keys); // sample X fragment
		float[][] transposedDataMatrix=General.transposeMatrix(dataMatrix); // fragment X sample
		float[] medianRelativeIntensities=new float[transposedDataMatrix.length];
		float[] intensityCVs=new float[medianRelativeIntensities.length];
		ArrayList<XYZPoint> points=new ArrayList<>(); 
		for (int i = 0; i < medianRelativeIntensities.length; i++) {
			medianRelativeIntensities[i]=QuickMedian.median(transposedDataMatrix[i]);
			intensityCVs[i]=General.stdev(transposedDataMatrix[i])/medianRelativeIntensities[i];
			points.add(new XYZPoint(-intensityCVs[i], medianRelativeIntensities[i], possibleMasses[i])); // negative CV (high is good)
		}
		Collections.sort(points);
		
		// return the best fragments reverse sorted on negative CV, intensity, then mass
		TDoubleArrayList d=new TDoubleArrayList();
		for (int i = points.size()-1; i >=0; i--) {
			d.add(points.get(i).z);
			if (d.size()>=targetNumberOfPoints) break;
			if (d.size()>=minNumberOfPoints&&points.get(i).x<-0.5f) break; // break at 50% CV if we've got enough peaks 
		}
		
		boolean[] usePeaks=new boolean[possibleMasses.length];
		for (int i = 0; i < usePeaks.length; i++) {
			if (d.contains(possibleMasses[i])) {
				usePeaks[i]=true;
			}
		}
		
		TObjectFloatHashMap<String> results=new TObjectFloatHashMap<>();
		for (int i = 0; i < dataMatrix.length; i++) {
			String sampleName=keys.get(i);
			float tic=0.0f;
			for (int j = 0; j < usePeaks.length; j++) {
				if (usePeaks[j]) {
					tic+=dataMatrix[i][j];
				}
			}
			results.put(sampleName, tic);
		}

		return d.toArray();
	}

	private Optional<double[]> checkEdgeCases(int n) {
		if (entryMap.size()==0) {
			return Optional.of(new double[0]);
		} else if (entryMap.size()==1) {
			// if there's only one point, then return the masses reverse sorted on correlation, then intensity, then mass
			LibraryEntry entry=entryMap.values().iterator().next();
			double[] masses=entry.getMassArray();
			float[] intensity=entry.getIntensityArray();
			float[] correlation=entry.getCorrelationArray();
			
			ArrayList<XYZPoint> points=new ArrayList<>();
			for (int i = 0; i < correlation.length; i++) {
				points.add(new XYZPoint(correlation[i], intensity[i], masses[i]));
			}
			Collections.sort(points);
			TDoubleArrayList d=new TDoubleArrayList();
			float sumIntensity=0.0f;
			for (int i = points.size()-1; i >=0; i--) {
				d.add(points.get(i).z);
				sumIntensity+=(float)points.get(i).y;
				if (d.size()>=n) break;
			}
			
			TObjectFloatHashMap<String> result=new TObjectFloatHashMap<>();
			result.put(entryMap.keySet().iterator().next(), sumIntensity);
			return Optional.of(d.toArray());
		}
		return Optional.empty();
	}

	private float[][] generateDataMatrix(double[] possibleMasses, ArrayList<String> keys) {
		float[][] dataMatrix=new float[keys.size()][]; // sample X fragment
		int count=0;
		for (String key : keys) {
			LibraryEntry entry=entryMap.get(key);
			
			if (entry!=null) {
				float[] integratedIntensities = getPossibleIntensities(possibleMasses, entry);
				dataMatrix[count]=integratedIntensities;
			} else {
				dataMatrix[count]=new float[possibleMasses.length]; // no intensities
			}
			count++;
		}
		return dataMatrix;
	}

	private float[] getPossibleIntensities(double[] possibleMasses, LibraryEntry entry) {
		double[] masses=entry.getMassArray();
		float[] intensities=entry.getIntensityArray();

		float[] integratedIntensities=tolerance.getIntegratedIntensities(masses, intensities, possibleMasses);
		integratedIntensities=General.normalize(integratedIntensities);
		return integratedIntensities;
	}
}
