package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import edu.washington.gs.maccoss.encyclopedia.utils.Quadruplet;
import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class PeakChromatogram extends Peak {
	private final float correlation;
	private final boolean isQuantified;

	public PeakChromatogram(double mass, float intensity, float correlation, boolean isQuantified) {
		super(mass, intensity);
		this.correlation=correlation;
		this.isQuantified=isQuantified;
	}
	
	public float getCorrelation() {
		return correlation;
	}
	
	public static Quadruplet<double[], float[], float[], boolean[]> toChromatogramArrays(Collection<PeakChromatogram> peaks) {
		double[] masses=new double[peaks.size()];
		float[] intensities=new float[peaks.size()];
		float[] correlations=new float[peaks.size()];
		boolean[] isQuantifiedArray=new boolean[peaks.size()];
		int count=0;
		for (PeakChromatogram peak : peaks) {
			masses[count]=peak.mass;
			intensities[count]=peak.intensity;
			correlations[count]=peak.correlation;
			isQuantifiedArray[count]=peak.isQuantified;
			count++;
		}
		return new Quadruplet<double[], float[], float[], boolean[]>(masses, intensities, correlations, isQuantifiedArray);
	}
	
	public static Quadruplet<double[], float[], float[], boolean[]> toChromatogramArrays(PeakChromatogram[] peaks) {
		double[] masses=new double[peaks.length];
	float[] intensities=new float[peaks.length];
	float[] correlations=new float[peaks.length];
	boolean[] isQuantifiedArray=new boolean[peaks.length];
	int count=0;
	for (PeakChromatogram peak : peaks) {
		masses[count]=peak.mass;
		intensities[count]=peak.intensity;
		correlations[count]=peak.correlation;
		isQuantifiedArray[count]=peak.isQuantified;
		count++;
	}
	return new Quadruplet<double[], float[], float[], boolean[]>(masses, intensities, correlations, isQuantifiedArray);
	}
	
	public static ArrayList<PeakChromatogram> fromChromatogramArrays(double[] masses, float[] intensities, float[] correlations, boolean[] isQuantifiedArray) {
		ArrayList<PeakChromatogram> peaks=new ArrayList<>();
		for (int i=0; i<masses.length; i++) {
			peaks.add(new PeakChromatogram(masses[i], intensities[i], correlations[i], isQuantifiedArray[i]));
		}
		return peaks;
	}
	
}
