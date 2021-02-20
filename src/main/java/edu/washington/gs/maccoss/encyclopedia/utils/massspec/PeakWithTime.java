package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.ArrayList;
import java.util.Collection;

import edu.washington.gs.maccoss.encyclopedia.utils.Triplet;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;

public class PeakWithTime extends Peak {
	private final float rtInSec;

	public PeakWithTime(double mass, float intensity, float rtInSec) {
		super(mass, intensity);
		this.rtInSec=rtInSec;
	}
	
	public float getRtInSec() {
		return rtInSec;
	}
	
	public static Triplet<double[], float[], float[]> toRTArrays(Collection<PeakWithTime> peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		TFloatArrayList rts=new TFloatArrayList();
		for (PeakWithTime peak : peaks) {
			masses.add(peak.mass);
			intensities.add(peak.intensity);
			rts.add(peak.rtInSec);
		}
		return new Triplet<double[], float[], float[]>(masses.toArray(), intensities.toArray(), rts.toArray());
	}
	
	public static Triplet<double[], float[], float[]> toRTArrays(PeakWithTime[] peaks) {
		TDoubleArrayList masses=new TDoubleArrayList();
		TFloatArrayList intensities=new TFloatArrayList();
		TFloatArrayList rts=new TFloatArrayList();
		for (PeakWithTime peak : peaks) {
			masses.add(peak.mass);
			intensities.add(peak.intensity);
			rts.add(peak.rtInSec);
		}
		return new Triplet<double[], float[], float[]>(masses.toArray(), intensities.toArray(), rts.toArray());
	}
	
	public static ArrayList<PeakWithTime> fromChromatogramArrays(double[] masses, float[] intensities, float[] rtsInSec) {
		ArrayList<PeakWithTime> peaks=new ArrayList<>();
		for (int i=0; i<masses.length; i++) {
			peaks.add(new PeakWithTime(masses[i], intensities[i], rtsInSec[i]));
		}
		return peaks;
	}

}
