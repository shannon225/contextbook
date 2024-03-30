package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Comparator;

public class PeakIntensityComparator implements Comparator<Peak> {
	public static final PeakIntensityComparator DEFAULT_INTENSITY_COMPARATOR=new PeakIntensityComparator();
	
	public PeakIntensityComparator() {
	}

	@Override
	public int compare(Peak o1, Peak o2) {
		if (o1==null&&o2==null) return 0;
		if (o1==null) return -1;
		if (o2==null) return 1;
		
		return Float.compare(o1.intensity, o2.intensity);
	}

}
