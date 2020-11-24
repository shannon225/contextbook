package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

import java.util.Comparator;

import edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation.TransitionRefiner;

public class QuantitativePeakIntensityComparator implements Comparator<PeakChromatogram> {
	public QuantitativePeakIntensityComparator() {
	}

	@Override
	public int compare(PeakChromatogram o1, PeakChromatogram o2) {
		if (o1==null&&o2==null) return 0;
		if (o1==null) return -1;
		if (o2==null) return 1;
		
		byte status1=getStatus(o1);
		byte status2=getStatus(o2);
		
		int c=Byte.compare(status1, status2);
		if (c!=0) return c;
		
		return Float.compare(o1.intensity, o2.intensity);
	}

	private byte getStatus(PeakChromatogram o1) {
		if (o1.getCorrelation()>TransitionRefiner.quantitativeCorrelationThreshold) {
			return 3;
		} else if (o1.getCorrelation()>TransitionRefiner.identificationCorrelationThreshold) {
			return 2;
		} else {
			return 1;
		}
	}

}
