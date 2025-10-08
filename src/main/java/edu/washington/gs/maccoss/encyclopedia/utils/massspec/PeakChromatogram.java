package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collection;

import edu.washington.gs.maccoss.encyclopedia.utils.Quadruplet;

public class PeakChromatogram extends Peak {
	private final float correlation;
	private final boolean isQuantified;

	public PeakChromatogram(double mass, float intensity, float correlation, boolean isQuantified) {
		super(mass, intensity);
		this.correlation=correlation;
		this.isQuantified=isQuantified;
	}
	
	public static PeakChromatogram mergePeaks(PeakChromatogram p1, PeakChromatogram p2) {
		float weightedSum=p1.intensity+p2.intensity;
		if (weightedSum<=0) weightedSum=1; // guard rail to not divide by 0
		double mass=(p1.mass*p1.intensity+p2.mass*p2.intensity)/weightedSum;
		return new PeakChromatogram(mass, weightedSum, Math.max(p1.correlation, p2.correlation), p1.isQuantified||p2.isQuantified);
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
