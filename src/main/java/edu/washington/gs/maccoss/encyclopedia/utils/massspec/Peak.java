package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Collection;

import edu.washington.gs.maccoss.encyclopedia.utils.Pair;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class Peak implements Comparable<Peak> {
	public final double mass;
	public final float intensity;
	
	public Peak(double mass, float intensity) {
		this.mass = mass;
		this.intensity = intensity;
	}
	
	public int compareTo(Peak o) {
		if (o==null) return 1;
		if (mass>o.mass) return 1;
		if (mass<o.mass) return -1;
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
}
