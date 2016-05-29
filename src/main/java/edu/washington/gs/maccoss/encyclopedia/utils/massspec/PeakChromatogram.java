package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Collection;

import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class PeakChromatogram extends Peak {
	private final float correlation;

	public PeakChromatogram(double mass, float intensity, float correlation) {
		super(mass, intensity);
		this.correlation=correlation;
	}
	
	public static Triplet<double[], float[], float[]> toChromatogramArrays(Collection<PeakChromatogram> peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		TFloatArrayList correlations=new TFloatArrayList();
		for (PeakChromatogram peak : peaks) {
			masses.add(peak.mass);
			intensities.add(peak.intensity);
			correlations.add(peak.correlation);
		}
		return new Triplet<double[], float[], float[]>(masses.toArray(), intensities.toArray(), correlations.toArray());
	}
	
	public static Triplet<double[], float[], float[]> toChromatogramArrays(PeakChromatogram[] peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		TFloatArrayList correlations=new TFloatArrayList();
		for (PeakChromatogram peak : peaks) {
			masses.add(peak.mass);
			intensities.add(peak.intensity);
			correlations.add(peak.correlation);
		}
		return new Triplet<double[], float[], float[]>(masses.toArray(), intensities.toArray(), correlations.toArray());
	}
	
}
