package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Collection;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import edu.washington.gs.maccoss.encyclopedia.utils.graphing.PointInterface;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class Peak implements PointInterface {
	public final double mass;
	public final float intensity;
	
	public Peak(double mass, float intensity) {
		this.mass = mass;
		this.intensity = intensity;
	}
	
	@Override
	public String toString() {
		return "("+mass+","+intensity+")";
	}

	@Override
	public double getX() {
		return mass;
	}

	@Override
	public double getY() {
		return intensity;
	}
	
	public float getIntensity() {
		return intensity;
	}
	public double getMass() {
		return mass;
	}
	
	/**
	 * doesn't compare on Y (intensity)
	 */
	@Override
	public int compareTo(PointInterface o) {
		if (o==null) return 1;
		if (mass>o.getX()) return 1;
		if (mass<o.getX()) return -1;
		return 0;
	}
	
	public static Pair<double[], float[]> toArrays(Collection<Peak> peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		for (Peak peak : peaks) {
			masses.add(peak.mass);
			intensities.add(peak.intensity);
		}
		return new Pair<double[], float[]>(masses.toArray(), intensities.toArray());
	}
	
	public static Pair<double[], float[]> toArrays(Peak[] peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		for (Peak peak : peaks) {
			masses.add(peak.mass);
			intensities.add(peak.intensity);
		}
		return new Pair<double[], float[]>(masses.toArray(), intensities.toArray());
	}
}
