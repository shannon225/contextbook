package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Comparator;

public class QuantitativePeakCorrelationComparator implements Comparator<PeakChromatogram> {
	public QuantitativePeakCorrelationComparator() {
	}

	@Override
	public int compare(PeakChromatogram o1, PeakChromatogram o2) {
		if (o1==null&&o2==null) return 0;
		if (o1==null) return -1;
		if (o2==null) return 1;
		
		int c=Float.compare(o1.getCorrelation(), o2.getCorrelation());
		if (c!=0) return c;
		
		c=Float.compare(o1.intensity, o2.intensity);
		if (c!=0) return c;
		
		return Double.compare(o1.mass, o2.mass);
	}

}
